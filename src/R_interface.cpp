#include "config.hpp"

#include <cstddef> // size_t
#include <cstdint>
#include <cstdlib> // malloc, free for PoDs
#include <cstring> // memcpy
#include <math.h> // nan


// R headers
#include <external/R.h>

#include <Rversion.h>

#if R_VERSION >= R_Version(3, 6, 2)
#define USE_FC_LEN_T
#endif

#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Random.h> // GetRNGstate, unif_rand

#undef USE_FC_LEN_T

// the flat C API of the installed dbarts (dbarts.h); the stubs resolve each
// entry point through R_GetCCallable on first use
#define DBARTS_USE_STUBS
#include <dbarts/dbarts.h>

#include "sensitivityAnalysis.hpp"
#include "guessNumCores.hpp"

#include "treatmentModel.hpp"
#include "probitTreatmentModel.hpp"
#include "probitEMTreatmentModel.hpp"
#include "bartTreatmentModel.hpp"

#include <glm/glm.h>

#include <external/random.h>

#include <misc/simd.h>

using std::size_t;
using std::uint32_t;
using std::strcmp;

namespace {
  SEXP getListElement(SEXP list, const char *str);
    
  enum TreatmentModelType {
    PROBIT_EM,
    PROBIT,
    BART
  };
  
  cibart::TreatmentModel* createTreatmentModel(SEXP modelExpr, TreatmentModelType* modelType,
                                               SEXP propControlExpr, SEXP propModelExpr, SEXP propDataExpr)
  {
    SEXP classExpr = GET_CLASS(modelExpr);
    if (isNull(classExpr) || !IS_CHARACTER(classExpr)) Rf_error("treatment model not of appropriate class");
    
    const char* className = CHAR(STRING_ELT(classExpr, 0));
    if (strcmp(className, "probitEMTreatmentModel") == 0) {
      *modelType = PROBIT_EM;
      return new cibart::ProbitEMTreatmentModel(static_cast<size_t>(INTEGER(getListElement(modelExpr, "maxIter"))[0]));
    }
    
    if (strcmp(className, "probitTreatmentModel") == 0) {
      *modelType = PROBIT;
      
      SEXP familyExpr = getListElement(modelExpr, "family");
      if (isNull(familyExpr) || !IS_CHARACTER(familyExpr)) Rf_error("probit treatment model lacks family attribute");
      
      cibart::ProbitPriorType family;
      const char* familyName = CHAR(STRING_ELT(familyExpr, 0));
      if (strcmp(familyName, "studentt") == 0) {
        family = cibart::PROBIT_PRIOR_STUDENT_T;
      } else if (strcmp(familyName, "normal") == 0) {
        family = cibart::PROBIT_PRIOR_NORMAL;
      } else if (strcmp(familyName, "flat") == 0) {
        family = cibart::PROBIT_PRIOR_FLAT;
      } else Rf_error("unrecognized probit treatment model family: '%s'", familyName);
      
      cibart::ProbitPrior* prior = NULL;
      switch(family) {
        case cibart::PROBIT_PRIOR_STUDENT_T:
        // prior = new cibart::ProbitStudentTPrior;
        // these are PoDs which require that we malloc them
        prior = static_cast<cibart::ProbitPrior*>(malloc(sizeof(cibart::ProbitStudentTPrior)));
        static_cast<cibart::ProbitStudentTPrior *>(prior)->scale = REAL(getListElement(modelExpr, "scale"));
        static_cast<cibart::ProbitStudentTPrior *>(prior)->dof   = REAL(getListElement(modelExpr, "df"))[0];
        break;
        case cibart::PROBIT_PRIOR_NORMAL:
        // prior = new cibart::ProbitNormalPrior;
        prior = static_cast<cibart::ProbitPrior*>(malloc(sizeof(cibart::ProbitNormalPrior)));
        static_cast<cibart::ProbitNormalPrior *>(prior)->scale = REAL(getListElement(modelExpr, "scale"));
        case cibart::PROBIT_PRIOR_FLAT:
        break;
      }
      
      return new cibart::ProbitTreatmentModel(family, prior);
    }
    
    if (strcmp(className, "bartTreatmentModel") == 0) {
      *modelType = BART;

      // the propensity model's dbarts spec triple is built R-side; ntree /
      // keepevery / k are baked into it, so nothing is read off modelExpr here
      if (propControlExpr == R_NilValue || propModelExpr == R_NilValue || propDataExpr == R_NilValue)
        Rf_error("bart treatment model requires its dbarts specification");

      return new cibart::BARTTreatmentModel(propControlExpr, propModelExpr, propDataExpr);
    }
    
    Rf_error("unrecognized treatment model: '%s'", className);
  }
  
  void destroyTreatmentModel(cibart::TreatmentModel* treatmentModelPtr, TreatmentModelType modelType)
  {
    switch(modelType) {
      case PROBIT_EM:
      delete static_cast<cibart::ProbitEMTreatmentModel*>(treatmentModelPtr);
      break;
      case PROBIT:
      {
        cibart::ProbitTreatmentModel* treatmentModel = static_cast<cibart::ProbitTreatmentModel*>(treatmentModelPtr);
        // delete treatmentModel->prior;
        free(const_cast<cibart::ProbitPrior*>(treatmentModel->prior));
        delete treatmentModel;
      }
      break;
      case BART:
      {
        delete static_cast<cibart::BARTTreatmentModel*>(treatmentModelPtr);
      }
      break;
    }
  }
  
  SEXP glmFit(SEXP y, SEXP n, SEXP x, SEXP w, SEXP offset)
  {
    double* yPtr = NULL;
    size_t numObs = static_cast<size_t>(XLENGTH(y));
    if (!isReal(y)) {
      if (!isInteger(y)) Rf_error("y must be of type real or integer");
      yPtr = new double[numObs];
      int* yInt = INTEGER(y);
      
      for (size_t i = 0; i < numObs; ++i) 
        yPtr[i] = (yInt[i] != 0.0 ? 1.0 : 0.0);
    }
      
    if (!isReal(x)) Rf_error("x must be of type real");
    
    int* dims;
    
    SEXP dimsExpr = GET_DIM(x);
    if (length(dimsExpr) != 2) Rf_error("x must be a matrix");
    dims = INTEGER(dimsExpr);
    if (static_cast<size_t>(dims[0]) != numObs) Rf_error("num rows in x must match length of y");
    
    size_t numCoefs = static_cast<size_t>(dims[1]);
    
    double* nPtr = NULL;
    if (!isNull(n)) {
      if (static_cast<size_t>(XLENGTH(n)) != numObs) Rf_error("length of n must match that of y");
      nPtr = REAL(n);
    }
    
    double* wPtr = NULL;
    if (!isNull(w)) {
      if (static_cast<size_t>(XLENGTH(w)) != numObs) Rf_error("length of w must match that of y");
      wPtr = REAL(w);
    }
    
    double* offsetPtr = NULL;
    if (!isNull(offset)) {
      if (static_cast<size_t>(XLENGTH(offset)) != numObs) Rf_error("length of offset must match that of y");
      offsetPtr = REAL(offset);
    }
    
    size_t scratchSize = glm_getDoubleScratchSize(numObs, numCoefs);
    double* scratch = new double[scratchSize];
    
    SEXP result = PROTECT(allocVector(REALSXP, static_cast<R_xlen_t>(numCoefs)));
    
    glm_fitGeneralizedLinearModel(yPtr == NULL ? REAL(y) : yPtr, nPtr, numObs, REAL(x), numCoefs, wPtr, offsetPtr,
                                  REAL(result), GLM_FAMILY_BINOMIAL, GLM_LINK_PROBIT, 30, scratch);
    
    delete [] scratch;
    if (yPtr != NULL) delete [] yPtr;
    
    UNPROTECT(1);
    
    return(result);
  }
  
  SEXP fitSensitivityAnalysis(SEXP y, SEXP z, SEXP xExpr,
                              SEXP x_testExpr,
                              SEXP zetaY, SEXP zetaZ, SEXP theta,
                              SEXP estimandExpr, SEXP treatmentModelExpr,
                              SEXP sensControl, SEXP verboseExpr,
                              SEXP outcomeControlExpr, SEXP outcomeModelExpr, SEXP outcomeDataExpr,
                              SEXP propControlExpr, SEXP propModelExpr, SEXP propDataExpr)
  {
    int* dims;
    
    if (!isReal(y)) Rf_error("y must be of type real.");
    if (!isReal(z)) Rf_error("z must be of type real.");
    if (!isReal(xExpr)) Rf_error("x must be of type real.");
    
    size_t numObservations = static_cast<size_t>(XLENGTH(y));
    if (numObservations == 0) Rf_error("y must have positive length");
    if (static_cast<size_t>(XLENGTH(z)) != numObservations) Rf_error("length of z and y must be equal");
    
    size_t numPredictors = 0;
    const double* x = NULL;
    if (XLENGTH(xExpr) > 0) {
      dims = INTEGER(getAttrib(xExpr, R_DimSymbol));
      if (dims == NULL || XLENGTH(getAttrib(xExpr, R_DimSymbol)) != 2) Rf_error("x must be a matrix");
      if (static_cast<size_t>(dims[0]) != numObservations) Rf_error("num rows of x and length of y must be equal");
      
      numPredictors = static_cast<size_t>(dims[1]);
      x = REAL(xExpr);
    }
    
    if (!isReal(x_testExpr)) Rf_error("x_test must be of type real.");
    
    size_t numTestObservations = 0;
    const double* x_test = NULL;
    if (XLENGTH(x_testExpr) > 0) {
      dims = INTEGER(getAttrib(x_testExpr, R_DimSymbol));
      if (dims == NULL || XLENGTH(getAttrib(x_testExpr, R_DimSymbol)) != 2) Rf_error("x_test must be a matrix");
      if (static_cast<size_t>(dims[1]) != numPredictors + 1) Rf_error("num columns of x_test must be equal to that of the column-combined z vector and x matrix");

      numTestObservations = static_cast<size_t>(dims[0]);
      x_test = REAL(x_testExpr);
    }
    
    if (!isReal(zetaY)) Rf_error("zetaY must be of type real");
    if (!isReal(zetaZ)) Rf_error("zetaZ must be of type real");
    
    size_t numZetaY = static_cast<size_t>(XLENGTH(zetaY));
    size_t numZetaZ = static_cast<size_t>(XLENGTH(zetaZ));
    
    if (!isString(estimandExpr)) Rf_error("estimand must be of type char");
    
    // turn estimand name into enum
    const char* estimandName = CHAR(STRING_ELT(estimandExpr, 0));
    cibart::EstimandType estimand;
    if (strncmp(estimandName, "ATE", 4) == 0) {
      estimand = cibart::ATE;
    } else if (strncmp(estimandName, "ATT", 4) == 0) {
      estimand = cibart::ATT;
    } else if (strncmp(estimandName, "ATC", 4) == 0) {
      estimand = cibart::ATC;
    } else {
      Rf_error("Illegal estimand type: %s. Must be 'ATE', 'ATT', or 'ATC'", estimandName);
    }
    
    TreatmentModelType treatmentModelType;
    cibart::TreatmentModel* treatmentModel = createTreatmentModel(treatmentModelExpr, &treatmentModelType,
                                                                  propControlExpr, propModelExpr, propDataExpr);
    
    SEXP sensParameterExpr = getListElement(sensControl, "n.sim");
    if (sensParameterExpr == R_NilValue) Rf_error("n.sim must be specified in iteration control");
    size_t numSimsPerCell = static_cast<size_t>(INTEGER(sensParameterExpr)[0]);
    
    sensParameterExpr = getListElement(sensControl, "n.burn.init");
    if (sensParameterExpr == R_NilValue) Rf_error("n.burn.init must be specified in iteration control");
    size_t numInitialBurnIn = static_cast<size_t>(INTEGER(sensParameterExpr)[0]);
    
    sensParameterExpr = getListElement(sensControl, "n.burn.cell");
    if (sensParameterExpr == R_NilValue) Rf_error("n.burn.cell must be specified in iteration control");
    size_t numCellSwitchBurnIn = static_cast<size_t>(INTEGER(sensParameterExpr)[0]);
    
    sensParameterExpr = getListElement(sensControl, "n.thin");
    if (sensParameterExpr == R_NilValue) Rf_error("n.thin must be specified in iteration control");
    size_t numTreeSamplesToThin = static_cast<size_t>(INTEGER(sensParameterExpr)[0]);
    
    // n.thread is honored no longer: the flat C API drives BART on the main R
    // thread only, so the grid runs sequentially (the classic C-level grid
    // parallelism relied on the removed per-thread RNG injection)

    if (!isLogical(verboseExpr)) Rf_error("verbose must be of type logical");
    if (length(verboseExpr) == 0) Rf_error("verbose must be of length at least 1");
    bool verbose = LOGICAL(verboseExpr)[0] != 0;
     
    
    // create result storage and make it user friendly
    SEXP dimsExpr, namesExpr;
    
    SEXP result = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(result, 0, allocVector(REALSXP, static_cast<R_xlen_t>(numZetaY * numZetaZ * numSimsPerCell)));
    SET_VECTOR_ELT(result, 1, allocVector(REALSXP, static_cast<R_xlen_t>(numZetaY * numZetaZ)));
    
    SEXP fittedCoefficients = VECTOR_ELT(result, 0);
    dimsExpr = PROTECT(dimsExpr = allocVector(INTSXP, 3));
    dims = INTEGER(dimsExpr);
    dims[0] = static_cast<int>(numSimsPerCell);
    dims[1] = static_cast<int>(numZetaY);
    dims[2] = static_cast<int>(numZetaZ);
    setAttrib(fittedCoefficients, R_DimSymbol, dimsExpr);

    
    SEXP standardErrors = VECTOR_ELT(result, 1);
    dimsExpr = PROTECT(dimsExpr = allocVector(INTSXP, 2));
    dims = INTEGER(dimsExpr);
    dims[0] = static_cast<int>(numZetaY);
    dims[1] = static_cast<int>(numZetaZ);
    setAttrib(standardErrors, R_DimSymbol, dimsExpr);
        
    
    setAttrib(result, R_NamesSymbol, namesExpr = allocVector(STRSXP, 2));
    SET_STRING_ELT(namesExpr, 0, mkChar("sens.coef"));
    SET_STRING_ELT(namesExpr, 1, mkChar("sens.se"));
    
    // draw the confounder-generator seed from R's stream in a narrow bracket;
    // the dbarts_sampler_* entry points manage R's RNG internally afterwards, so
    // no outer bracket may span them (it would double-bracket R's stream)
    GetRNGstate();
    uint_least32_t rngSeed = static_cast<uint_least32_t>(unif_rand() * 4294967295.0);
    PutRNGstate();

    cibart::fitSensitivityAnalysis(REAL(y), REAL(z), x,
                                   numObservations, numPredictors,
                                   x_test,
                                   numTestObservations,
                                   REAL(zetaY), REAL(zetaZ), numZetaY, numZetaZ,
                                   REAL(theta)[0],
                                   estimand, *treatmentModel,
                                   numSimsPerCell,
                                   numInitialBurnIn,
                                   numCellSwitchBurnIn,
                                   numTreeSamplesToThin,
                                   outcomeControlExpr, outcomeModelExpr, outcomeDataExpr,
                                   rngSeed,
                                   REAL(fittedCoefficients),
                                   REAL(standardErrors),
                                   verbose);

    UNPROTECT(3);
    
    destroyTreatmentModel(treatmentModel, treatmentModelType);
    
    return result;
  }
  
  SEXP getListElement(SEXP list, const char* str)
  {
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    for (int i = 0; i < length(list); i++)
      if (std::strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
        elmt = VECTOR_ELT(list, i);
        break;
      }
    return elmt;
  }
  
  SEXP guessNumCores()
  {
    uint32_t numPhyiscalProcessors, numLogicalProcessors;
    cibart::guessNumCores(&numPhyiscalProcessors, &numLogicalProcessors);
    
    SEXP resultExpr = allocVector(INTSXP, 2);
    PROTECT(resultExpr);
    int* result = INTEGER(resultExpr);
    
    result[0] = numPhyiscalProcessors <= 0 ? NA_INTEGER : static_cast<int>(numPhyiscalProcessors);
    result[1] = numLogicalProcessors  <= 0 ? NA_INTEGER : static_cast<int>(numLogicalProcessors);
    UNPROTECT(1);
    
    return resultExpr;
  }
}

#if __cplusplus >= 202002L
#  include <bit>
#else

#  if __cplusplus >= 201103L
#    include <type_traits>

namespace std {

template <class To, class From>
typename std::enable_if<
  sizeof(To) == sizeof(From) &&
  std::is_trivially_copyable<From>::value &&
  std::is_trivially_copyable<To>::value,
  To>::type
// constexpr support needs compiler magic
bit_cast(const From& src) noexcept
{
  static_assert(std::is_trivially_constructible<To>::value,
    "This implementation additionally requires destination type to be trivially constructible");

  To dst;
  std::memcpy(&dst, &src, sizeof(To));
  return dst;
}

#  else

// We are only using this to cast function pointers, which are trivially copiable.
// is_trivially_copyable is compiler specific and isn't worth trying to drop in
// an implementation.
template <class To, class From>
To
bit_cast(const From& src)
{
  To dst;
  std::memcpy(&dst, &src, sizeof(To));
  return dst;
}

#  endif

}

#endif

extern "C" {

#define DEF_FUNC(_N_, _F_, _A_) { _N_, std::bit_cast<DL_FUNC>(&_F_), _A_ }
  
  R_CallMethodDef R_callMethods[] = {
    DEF_FUNC("treatSens_fitSensitivityAnalysis", fitSensitivityAnalysis, 17),
    DEF_FUNC("treatSens_guessNumCores", guessNumCores, 0),
    DEF_FUNC("treatSens_glmFit", glmFit, 5),
    {NULL, NULL, 0}
  };
}

#undef DEF_FUNC

extern "C" {
  
  void R_init_treatSens(DllInfo* info)
  {
    R_registerRoutines(info, NULL, R_callMethods, NULL, NULL);
    R_useDynamicSymbols(info, static_cast<Rboolean>(FALSE));

    // dbarts flat C API handshake: major must match, minor must be at least ours
    //
    // Pre-release lockstep check: the version constants do not move before
    // 1.0-0, so the major/minor handshake alone can never fire and only the
    // exact signature token catches a stale binary. Remove or downgrade it at
    // the 1.0-0 freeze - post-release a legitimate minor append moves the hash
    // and a hard equality would refuse every consumer until it is rebuilt.
    if (dbarts_apiMajorVersion() != DBARTS_C_API_MAJOR || dbarts_apiMinorVersion() < DBARTS_C_API_MINOR ||
        dbarts_apiHash() != DBARTS_C_API_HASH)
      Rf_error("treatSens was built against dbarts C API %d.%d but the installed dbarts provides %d.%d, or the two carry different entry-point signatures; reinstall treatSens",
               DBARTS_C_API_MAJOR, DBARTS_C_API_MINOR, dbarts_apiMajorVersion(), dbarts_apiMinorVersion());

    misc_simd_init();
  }
}

