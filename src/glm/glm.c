#include <glm/glm.h>

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <float.h>

#include <external/io.h>
#include <external/stats.h>
#include <external/linearAlgebra.h>

#define stdNormalPDF(_X_) ext_densityOfNormal((_X_), 0.0, 1.0)
#define stdNormalCDF(_Q_) ext_cumulativeProbabilityOfNormal((_Q_), 0.0, 1.0)
#define stdNormalInvCDF(_P_) ext_quantileOfNormal((_P_), 0.0, 1.0)

#define LEAST_SQUARES_TOLERANCE 1.0e-7
#define DEVIANCE_TOLERANCE 1.0e-8

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

// function pointers
// takes a vector of doubles, applies func to all and stores result in 2nd arg
typedef void (*glm_doubleMapFunction)(const double* restrict, double* restrict, size_t);
// takes a vector of doubles, returns a bool
typedef bool (*glm_boolFunction)(const double* restrict, size_t);
// takes vectors of doubles y, E[y], and weights; returns a double
typedef double (*glm_doubleFunction)(const double* restrict, const double* restrict, const double*, size_t);

// used for storing said function pointers
typedef struct {
  glm_doubleMapFunction applyLink;
  glm_doubleMapFunction applyInverseLink;
  glm_doubleMapFunction applyDmuDeta; // point-wise application of derivative of mean function w/r/t linear predictor
  glm_doubleMapFunction applyVariance; // for model and link, calculates the variance as a func of mean
  
  glm_boolFunction isValidMean;
  glm_boolFunction isValidLinearPredictor;
  
  glm_doubleFunction getDevianceOfResiduals;
} glm_functions;


// forward declarations
static void initialize(double* y, double* weights, double* offset, double* mu,
                       const double* _y, const double* n, const double* _weights,
                       const double* _offset, size_t numObs);
static bool allFinite(const double* restrict array, size_t n);
static bool allPositive(const double* restrict array, size_t n);

static void getFunctionsFromFamilyAndLink(glm_family_t family, glm_link_t link, glm_functions* functions);

static void calculateZ(const double* eta, const double* offset, const double* y, const double* mu,
                       const double* dmudeta, const bool* goodObs, size_t numObs, double* z);
static void calculateObservationWeights(const double* priorWeights, const double* dmudeta,
                                        const double* variances, const bool* goodObs, size_t numObs,
                                        double* obsrvWeights);
static void subsetX(const double* _x, const bool* goodObs, size_t numObs, size_t numCoefs,
                    double* x, size_t numGoodObs);

static void weightZ(const double* restrict obsrvWeights, size_t numGoodObs, double* restrict z);
static void weightX(const double* restrict obsrvWeights, size_t numGoodObs, size_t numCoefs, double* restrict x);

size_t glm_getDoubleScratchSize(size_t numObs, size_t numCoefs)
{
  return numObs * (9 + numCoefs);
}

void glm_fitGeneralizedLinearModel(const double* _y, const double* n, size_t numObs,
                                   const double *_x, size_t numCoefs,
                                   const double* _priorWeights, const double* _offset,
                                   double* prevCoefs,
                                   glm_family_t family, glm_link_t link,
                                   uint32_t maxIterations, double* scratch)
{
  size_t index;
  uint32_t optimIter, walkbackIter;
  
  double* heapOHeap = scratch;
  if (heapOHeap == NULL) heapOHeap = malloc(glm_getDoubleScratchSize(numObs, numCoefs) * sizeof(double));
  if (heapOHeap == NULL) ext_throwError("Unable to allocate memory for glm fit.");
  
  // quite likely some memory could be reused
  double* mu      = heapOHeap;
  double* eta     = heapOHeap +     numObs;
  double* dmudeta = heapOHeap + 2 * numObs;
  
  double* y            = heapOHeap + 3 * numObs; // y / n, mostly
  double* priorWeights = heapOHeap + 4 * numObs; // fill in 1s, in case of NULL
  double* obsrvWeights = heapOHeap + 5 * numObs; // mus and eta everywhere
  double* offset       = heapOHeap + 6 * numObs;
  double* variances    = heapOHeap + 7 * numObs;
  
  double* z            = heapOHeap + 8 * numObs; // among other things, these may have dropped rows
  double* x            = heapOHeap + 9 * numObs; // I wish I had a better name for z, but it pops out
                                                 // when you do the glm math
  

  const char* errorMessage = NULL;
  double currCoefs[numCoefs];
  
  glm_functions functions;
  
  getFunctionsFromFamilyAndLink(family, link, &functions);
    
  glm_doubleMapFunction applyLink = functions.applyLink;
  glm_doubleMapFunction applyInverseLink = functions.applyInverseLink;
  glm_doubleMapFunction applyDmuDeta = functions.applyDmuDeta;
  glm_doubleMapFunction applyVariance = functions.applyVariance;
  
  glm_boolFunction isValidMean = functions.isValidMean;
  glm_boolFunction isValidLinearPredictor = functions.isValidLinearPredictor;
  
  glm_doubleFunction getDevianceOfResiduals = functions.getDevianceOfResiduals;
  
  
  uint32_t numGoodObs;
  bool goodObs[numObs];
  double currDeviance, prevDeviance;
  int32_t leastSquaresResult;
  char* leastSquaresMessage = NULL;
  
  // end of declarations
  
  
  initialize(y, priorWeights, offset, mu, _y, n, _priorWeights, _offset, numObs);
  
  applyLink((const double*) mu, eta, numObs);
  if (!isValidMean(mu, numObs) || !isValidLinearPredictor(eta, numObs)) {
    errorMessage = "Cannot find valid starting values."; goto glm_cleanup;
  }
    
  prevDeviance = getDevianceOfResiduals(y, mu, priorWeights, numObs);
  for (index = 0; index < numCoefs; ++index) prevCoefs[index] = 0.0;
  
  for (optimIter = 0; optimIter < maxIterations; ++optimIter) {
    applyVariance((const double*) mu, variances, numObs);
    
    if (!allFinite  (variances, numObs)) { errorMessage = "NAN or Inf variance generated.";   goto glm_cleanup; }
    if (!allPositive(variances, numObs)) { errorMessage = "Non-positive variance generated."; goto glm_cleanup; }
    
    
    applyDmuDeta((const double*) eta, dmudeta, numObs);
    
    if (!allFinite(dmudeta, numObs)) { errorMessage = "NAN or Inf in dmu/deta."; goto glm_cleanup; }
    
    
    numGoodObs = 0;
    for (index = 0; index < numObs; ++index) {
      goodObs[index] = (dmudeta[index] != 0.0);
      if (goodObs[index] == true) ++numGoodObs;
    }
    if (numGoodObs == 0) { errorMessage = "No observations informative."; goto glm_cleanup; }
    
    calculateZ((const double*) eta, (const double*) offset, (const double*) y, (const double*) mu,
               (const double*) dmudeta, (const bool*) goodObs, numObs, z);
    calculateObservationWeights((const double*) priorWeights, (const double*) dmudeta,
                                (const double*) variances, (const bool*) goodObs, numObs, obsrvWeights);
    
    subsetX((const double*) _x, (const bool*) goodObs, numObs, numCoefs,
            x, numGoodObs);
    
    weightZ((const double*) obsrvWeights, numGoodObs, z);
    weightX((const double*) obsrvWeights, numGoodObs, numCoefs, x);
    
    leastSquaresResult = ext_findLeastSquaresFit(z, numObs, x, numCoefs, currCoefs, LEAST_SQUARES_TOLERANCE, NULL, &leastSquaresMessage);
    if (leastSquaresResult == 0) {
      ext_printMessage("least squares failed with error: '%s' at iteration %d.", leastSquaresMessage, optimIter);
      break;
    } else if (leastSquaresResult < 0) {
      errorMessage = leastSquaresMessage;
      goto glm_cleanup;
    } else if ((int32_t) numObs < leastSquaresResult) {
      errorMessage = "X matrix of greater rank than num observations."; goto glm_cleanup;
    }
    
    // calulate linPred = X %*% beta + offset
    ext_leftMultiplyMatrixAndVector((const double*) _x, numObs, numCoefs, (const double*) currCoefs, eta);
    ext_addVectorsInPlace((const double*) offset, numObs, 1.0, eta);
    
    applyInverseLink((const double*) eta, mu, numObs);
    currDeviance = getDevianceOfResiduals(y, mu, priorWeights, numObs);
    
    
    bool coefsValid = !isnan(currDeviance) && !isinf(currDeviance) && isValidMean(mu, numObs) && isValidLinearPredictor(eta, numObs);
    for (walkbackIter = 0; !coefsValid && walkbackIter < maxIterations; ++walkbackIter) {
      // step half way back to previous iteration
      for (size_t index = 0; index < numCoefs; ++index) currCoefs[index] = 0.5 * (currCoefs[index] + prevCoefs[index]);
      
      ext_leftMultiplyMatrixAndVector((const double*) _x, numObs, numCoefs, (const double*) currCoefs, eta);
      for (size_t index = 0; index < numObs; ++index) eta[index] += offset[index];
      
      applyInverseLink((const double*) eta, mu, numObs);
      currDeviance = getDevianceOfResiduals(y, mu, priorWeights, numObs);
    }
    if (walkbackIter >= maxIterations) { errorMessage = "Cannot correct step size to find valid coefficients."; goto glm_cleanup; }
    
    
    for (index = 0; index < numCoefs; ++index) prevCoefs[index] = currCoefs[index];
    
    if (fabs(currDeviance - prevDeviance) / (0.1 + fabs(currDeviance)) < DEVIANCE_TOLERANCE)
      break;
    
    prevDeviance = currDeviance;
  }
    
glm_cleanup:
  if (scratch == NULL) free(heapOHeap);
  
  if (errorMessage != NULL) ext_throwError(errorMessage);
}

static void initialize(double* y, double* weights, double* offset, double* mu,
                       const double* inputY, const double* n, const double* inputWeights,
                       const double* inputOffset, size_t numObs)
{
  double n_i;
  for (size_t i = 0; i < numObs; ++i) {
    n_i = (n == NULL ? 1.0 : n[i]);
    y[i] = (n_i == 0.0 ? 0.0 : inputY[i] / n_i);
    weights[i] = (inputWeights != NULL ? inputWeights[i] * n_i : n_i);
    mu[i] = (n_i * y[i] + 0.5) / (n_i + 1.0);
    offset[i] = (inputOffset != NULL ? inputOffset[i] : 0.0);
  }
}

static bool allFinite(const double* restrict array, size_t n)
{
  for (size_t i = 0; i < n; ++i)
    if (isnan(array[i]) || isinf(array[i])) return false;
  
  return true;
}

// assumes are valid, finite numbers at this point
static bool allPositive(const double* restrict array, size_t n)
{
  for (size_t i = 0; i < n; ++i)
    if (array[i] <= 0.0) return false;
  
  return true;
}

static void calculateZ(const double* eta, const double* offset, const double* y, const double* mu,
                       const double* dmudeta, const bool* goodObs, size_t numObs, double* z)
{
  size_t targetRow = 0;
  for (size_t sourceRow = 0; sourceRow < numObs; ++sourceRow) {
    if (goodObs[sourceRow] == true) {
      z[targetRow++] = (eta[sourceRow] - offset[sourceRow]) +
                       (y[sourceRow] - mu[sourceRow]) / dmudeta[sourceRow];
    }
  }
  for (; targetRow < numObs; ++targetRow) z[targetRow] = 0.0;
}

static void calculateObservationWeights(const double* priorWeights, const double* dmudeta,
                                        const double* variances, const bool* goodObs, size_t numObs,
                                        double* obsrvWeights)
{
  size_t targetRow = 0;
  for (size_t sourceRow = 0; sourceRow < numObs; ++sourceRow) {
    if (goodObs[sourceRow] == true) {
      obsrvWeights[targetRow++] =
        sqrt((priorWeights[sourceRow] * dmudeta[sourceRow] * dmudeta[sourceRow]) / variances[sourceRow]);
    }
  }
  for (; targetRow < numObs; ++targetRow) obsrvWeights[targetRow] = 0.0;
}

static void subsetX(const double* _x, const bool* goodObs, size_t numObs, size_t numCoefs,
                    double* x, size_t numGoodObs)
{
  size_t targetRow = 0;
  // column-major storage makes subsetting rows a pain
  for (size_t sourceRow = 0; sourceRow < numObs; ++sourceRow) {
    if (!goodObs[sourceRow]) continue;
    
    for (size_t col = 0; col < numCoefs; ++col) {
      x[targetRow + col * numGoodObs] = _x[sourceRow + col * numObs];
    }
    ++targetRow;
  }
  
  size_t baseOffset = numGoodObs * numCoefs;
  size_t zeroFillEnd = (numObs  - numGoodObs) * numCoefs + baseOffset;
  for (size_t i = baseOffset; i < zeroFillEnd; ++i) x[i] = 0.0;
}

static void weightZ(const double* restrict obsrvWeights, size_t numObs, double* restrict z)
{
  for (size_t i = 0; i < numObs; ++i) z[i] *= obsrvWeights[i];
}
// rows of X are reweighted; X is stored column major
static void weightX(const double* restrict obsrvWeights, size_t numObs, size_t numCoefs, double* restrict x)
{
  size_t offset = 0;
  for (size_t col = 0; col < numCoefs; ++col) {
    for (size_t row = 0; row < numObs; ++row) {
      x[offset++] *= obsrvWeights[row];
    }
  }
}

// ----------------------
// definitions for probit
static void probitLink(const double* restrict source, double* restrict target, size_t n) {
  for (size_t i = 0; i < n; ++i) target[i] = stdNormalInvCDF(source[i]);
}

static void probitInvLink(const double* restrict source, double* restrict target, size_t n) {
  for (size_t i = 0; i < n; ++i) target[i] = stdNormalCDF(source[i]);
}

static void probitDmuDeta(const double* restrict source, double* restrict target, size_t n) {
  double temp;
  for (size_t i = 0; i < n; ++i) {
    temp = stdNormalPDF(source[i]);
    target[i] = (temp < DBL_EPSILON ? DBL_EPSILON : temp);
  }
}

static void binomialVariance(const double* restrict mean, double* restrict var, size_t n) {
  for (size_t i = 0; i < n; ++i) var[i] = mean[i] * (1.0 - mean[i]);
}

static bool probitValidMu(const double* restrict mu, size_t n) {
  for (size_t i = 0; i < n; ++i)
    if (mu[i] <= 0.0 || mu[i] >= 1.0) return false;

  return true;
}

static bool probitValidEta(UNUSED const double* restrict eta, UNUSED size_t n) {
  return true;
}

#define y_log_y(_Y_, _MU_) ((_Y_) != 0.0 ? (_Y_) * log((_Y_) / (_MU_)) : 0.0)

static double probitDevResid(const double* restrict y, const double* restrict mu,
                             const double* restrict weights, size_t n) {
  double result = 0.0;
  for (size_t i = 0; i < n; ++i) {
    result += 2.0 * weights[i] * (y_log_y(y[i], mu[i]) + y_log_y(1.0 - y[i], 1.0 - mu[i]));
  }
  return result;
}

static void getFunctionsFromFamilyAndLink(glm_family_t family, glm_link_t link, glm_functions* functions)
{
  if (family == GLM_FAMILY_BINOMIAL) {
    if (link == GLM_LINK_PROBIT) {
      functions->applyLink = &probitLink;
      functions->applyInverseLink = &probitInvLink;
      functions->applyDmuDeta = &probitDmuDeta;
      functions->applyVariance = &binomialVariance;
      
      functions->isValidMean = &probitValidMu;
      functions->isValidLinearPredictor = &probitValidEta;
      
      functions->getDevianceOfResiduals = &probitDevResid;
    } else {
      ext_throwError("only binomial, probit implemented at present");
    }
  } else ext_throwError("only binomial, probit implemented at present");
}
