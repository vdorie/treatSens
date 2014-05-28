#include "config.hpp"

#include <cstring>
#include <cstddef>
#include <bart/cstdint>

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include "sensitivityAnalysis.hpp"

#include <bart/bartFit.hpp>
#include <bart/control.hpp>
#include <bart/data.hpp>
#include <bart/model.hpp>
#include <bart/results.hpp>

using std::size_t;

extern "C" {
  static SEXP getListElement(SEXP list, const char *str);
  
  SEXP cibart_fitSensitivityAnalysis(SEXP y, SEXP z, SEXP xExpr,
                                     SEXP z_testExpr, SEXP x_testExpr,
                                     SEXP zetaY, SEXP zetaZ, SEXP theta,
                                     SEXP estimandExpr, SEXP confounderModelExpr,
                                     SEXP fitControl, SEXP verboseExpr)
  {
    int* dims;
    
    if (!isReal(y)) error("y must be of type real.");
    if (!isReal(z)) error("z must be of type real.");
    if (!isReal(xExpr)) error("x must be of type real.");
    
    size_t numObservations = length(y);
    if (numObservations == 0) error("y must have positive length.");
    if (length(z) != (int) numObservations) error("length of z and y must be equal.");
    
    size_t numPredictors = 0;
    const double* x = NULL;
    if (length(xExpr) > 0) {
      dims = INTEGER(getAttrib(xExpr, R_DimSymbol));
      if (dims == NULL || length(getAttrib(xExpr, R_DimSymbol)) != 2) error("x must be a matrix.");
      if (dims[0] != (int) numObservations) error("num rows of x and length of y must be equal.");
      
      numPredictors = dims[1];
      x = REAL(xExpr);
    }
    
    if (!isLogical(z_testExpr)) error("z_test must be of type logical.");
    if (!isReal(x_testExpr)) error("x_test must be of type real.");
    
    if (length(z_testExpr) != (int) numObservations) error("length of z_test and y must be equal.");
    bool* z_test = new bool[numObservations];
    int* intTemp = LOGICAL(z_testExpr);
    for (size_t i = 0; i < numObservations; ++i) z_test[i] = intTemp[i] != 0;
    
    size_t numTestObservations = 0;
    const double* x_test = NULL;
    if (length(x_testExpr) > 0) {
      dims = INTEGER(getAttrib(x_testExpr, R_DimSymbol));
      if (dims == NULL || length(getAttrib(x_testExpr, R_DimSymbol)) != 2) error("x_test must be a matrix.");
      if (dims[1] != (int) (numPredictors + 1)) error("num columns of x_test must be equal to that of the column-combined z vector and x matrix.");

      numTestObservations = dims[0];
      x_test = REAL(x_testExpr);
    }
    
    if (!isReal(zetaY)) error("zetaY must be of type real.");
    if (!isReal(zetaZ)) error("zetaZ must be of type real.");
    
    size_t numZetaY = length(zetaY);
    size_t numZetaZ = length(zetaZ);
    
    if (!isString(estimandExpr)) error("estimand must be of type char.");
    
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
      error("Illegal estimand type: %s. Must be 'ATE', 'ATT', or 'ATC'.", estimandName);
    }
    
    
    if (!isString(confounderModelExpr)) error("estimand must be of type char.");
    
    // turn model name into enum
    const char* modelName = CHAR(STRING_ELT(confounderModelExpr, 0));
    cibart::ConfounderModelType confounderModel;
    if (strncmp(modelName, "normal", 7) == 0) {
      confounderModel = cibart::NORMAL;
    } else if (strncmp(modelName, "binomial", 9) == 0) {
      confounderModel = cibart::BINOMIAL;
    } else {
      error("Illegal confounder model type: %s. Must be 'normal' or 'binomial'.", modelName);
    }
    
    if (confounderModel == cibart::NORMAL)
      error("BART sensitivity analysis only implemented for binary confounders.");
    
    SEXP fitParameterExpr = getListElement(fitControl, "nsim");
    if (fitParameterExpr == R_NilValue) error("nsim must be specified in iteration control.");
    size_t numSimsPerCell = INTEGER(fitParameterExpr)[0];

    fitParameterExpr = getListElement(fitControl, "nchain");
    if (fitParameterExpr == R_NilValue) error("nchain must be specified in iteration control.");
    size_t numChains = INTEGER(fitParameterExpr)[0];
    
    fitParameterExpr = getListElement(fitControl, "nburn");
    if (fitParameterExpr == R_NilValue) error("nburn must be specified in iteration control.");
    size_t numBurnIn = INTEGER(fitParameterExpr)[0];
    
    fitParameterExpr = getListElement(fitControl, "nthin");
    if (fitParameterExpr == R_NilValue) error("nthin must be specified in iteration control.");
    size_t numTreeSamplesToThin = INTEGER(fitParameterExpr)[0];
    
    fitParameterExpr = getListElement(fitControl, "nthread");
    if (fitParameterExpr == R_NilValue) error("nthread must be specified in iteration control.");
    size_t numThreads = INTEGER(fitParameterExpr)[0];
    
    if (!isLogical(verboseExpr)) error("verbose must be of type logical.");
    if (length(verboseExpr) == 0) error("verbose must be of length at least 1.");
    bool verbose = LOGICAL(verboseExpr)[0] != 0;
     
    
    // create result storage and make it user friendly
    SEXP dimsExpr, namesExpr;
    
    SEXP result = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(result, 0, allocVector(REALSXP, numZetaY * numZetaZ * numSimsPerCell));
    SET_VECTOR_ELT(result, 1, allocVector(REALSXP, numZetaY * numZetaZ));
    
    SEXP fittedCoefficients = VECTOR_ELT(result, 0);
    dimsExpr = PROTECT(dimsExpr = allocVector(INTSXP, 3));
    dims = INTEGER(dimsExpr);
    dims[0] = numSimsPerCell;
    dims[1] = numZetaY;
    dims[2] = numZetaZ;
    setAttrib(fittedCoefficients, R_DimSymbol, dimsExpr);

    
    SEXP standardErrors = VECTOR_ELT(result, 1);
    dimsExpr = PROTECT(dimsExpr = allocVector(INTSXP, 2));
    dims = INTEGER(dimsExpr);
    dims[0] = numZetaY;
    dims[1] = numZetaZ;
    setAttrib(standardErrors, R_DimSymbol, dimsExpr);
        
    
    setAttrib(result, R_NamesSymbol, namesExpr = allocVector(STRSXP, 2));
    SET_STRING_ELT(namesExpr, 0, mkChar("sens.coef"));
    SET_STRING_ELT(namesExpr, 1, mkChar("sens.se"));
    
    GetRNGstate();
    
    cibart::fitSensitivityAnalysis(REAL(y), REAL(z), x,
                                 numObservations, numPredictors,
                                 z_test, x_test,
                                 numTestObservations,
                                 REAL(zetaY), REAL(zetaZ), numZetaY, numZetaZ,
                                 REAL(theta)[0],
                                 estimand, confounderModel,
                                 numSimsPerCell,
                                 numChains,
                                 numBurnIn,
                                 numTreeSamplesToThin,
                                 numThreads,
                                 REAL(fittedCoefficients),
                                 REAL(standardErrors),
                                 verbose);
    
    PutRNGstate();
    
    UNPROTECT(3);
        
    delete [] z_test;
    
    return(result);
  }
  
  static SEXP getListElement(SEXP list, const char* str)
  {
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    for (int i = 0; i < length(list); i++)
      if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
        elmt = VECTOR_ELT(list, i);
        break;
      }
    return elmt;
  }
  
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
  
  static R_CallMethodDef callMethods[] = {
    CALLDEF(cibart_fitSensitivityAnalysis, 12),
    {NULL, NULL, 0}
  };
  
  void R_init_SensBART(DllInfo* info)
  {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  }
}
