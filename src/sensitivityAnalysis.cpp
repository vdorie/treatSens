#include "config.hpp"
#include "sensitivityAnalysis.hpp"

#include <cstring>
#include <cmath>
#include <dbarts/cstdint.hpp>
#if !defined(HAVE_SYS_TIME_H) && defined(HAVE_GETTIMEOFDAY)
#undef HAVE_GETTIMEOFDAY
#endif
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h> // gettimeofday
#else
#include <time.h>
#endif

#include <dbarts/bartFit.hpp>
#include <dbarts/types.hpp>
#include <dbarts/results.hpp>
#include <dbarts/R_C_interface.hpp>

#include <external/io.h>
#include <external/stats.h>
#include <external/linearAlgebra.h>

#include <glm/glm.h>

#include <external/alloca.h>

#include <R_ext/Rdynload.h>

#define DEFAULT_BART_MAX_NUM_CUTS 100u

using std::size_t;
using std::uint32_t;

namespace {
  using namespace cibart;
  struct Control {
    EstimandType estimand;
    ConfounderModelType confounderModel;
    
    size_t numSimsPerCell;
    size_t numInitialBurnIn;
    size_t numCellSwitchBurnIn;
    size_t numTreeSamplesToThin;
    size_t numThreads;
    
    double theta; // prior probability of U == 1
    
    bool verbose;
    
    void (*initializeFit)(dbarts::BARTFit* fit, dbarts::Control* control, dbarts::Model* model, dbarts::Data* data);
    void (*invalidateFit)(dbarts::BARTFit* fit);
    dbarts::Results* (*runSampler)(dbarts::BARTFit* fit);
    dbarts::Results* (*runSamplerForIterations)(dbarts::BARTFit* fit, size_t numBurnIn, size_t numSamples);
    void (*setResponse)(dbarts::BARTFit* fit, const double* newResponse);
    void (*initializeCGMPrior)(dbarts::CGMPrior* prior, double, double);
    void (*invalidateCGMPrior)(dbarts::CGMPrior* prior);
    void (*initializeNormalPrior)(dbarts::NormalPrior* prior, const dbarts::Control* control, double);
    void (*invalidateNormalPrior)(dbarts::NormalPrior* prior);
    void (*initializeChiSquaredPrior)(dbarts::ChiSquaredPrior* prior, double, double);
    void (*invalidateChiSquaredPrior)(dbarts::ChiSquaredPrior* prior);
  };
  
  struct Data {
    // stuff passed in
    const double* y;
    const double* z;
    const double* x;
    size_t numObservations;
    size_t numPredictors;
    
    const bool* z_test;
    const double* x_test;
    size_t numTestObservations;
        
    // various transformations
    const double* bart_x_train;
    size_t numBartPredictors;
    const double* glm_x;
    size_t numGLMPredictors;
    
    Data(const double* y, const double* z, const double* x,
         size_t numObservations, size_t numPredictors, const bool* z_test, const double* x_test,
         size_t numTestObservations);
    ~Data();
  };
  
  struct Scratch {
    double* yMinusZetaU;
    double* p;
    double* u;
    double* glmScratch;
    
    double* temp_numObs_1;
    
    Scratch(const Control& control, const Data& data);
    ~Scratch();
  };
  
  // forward declarations
  //void getEstimatesForGridPoint(const Control& control, const Data& data, Scratch& scratch, double zetaY, double zetaZ,
  //                              double* estimates);
  void sampleConfounders(const Data& data, Scratch& scratch);
  void subtractConfounderFromResponse(const Data& data, Scratch& scratch, double zetaY);
  double estimateSigma(const Data& data, Scratch& scratch);
  void updateConfounderProbabilities(const Control& control, const Data& data, Scratch& scratch, double zetaY, double zetaZ,
                                     const double* yMinusZetaUHat, double sigma);
  void estimateTreatmentEffect(const Control& control, const Data& data,
                               const double* trainingSamples, const double* testSamples, double* estimates);
  void lookupBARTFunctions(Control& control);
  void sensitivityAnalysisTask(void* v_data);
  
#ifdef HAVE_GETTIMEOFDAY
  double subtractTimes(struct timeval end, struct timeval start);
#else
  double subtractTimes(time_t end, time_t start);
#endif
  
  struct GridCell {
    double zetaY;
    double zetaZ;
    size_t offset;
    size_t cellNumber;
  };
  
  struct ThreadData {
    const GridCell* gridCells;
    size_t numGridCells;
    size_t totalNumGridCells;
    
    Control* control;
    Data* data;
    Scratch* scratch;
    
    double* estimates;
    double* standardErrors;
  };
}

namespace cibart {
  void
  fitSensitivityAnalysis(const double* y,     // numObs x 1
                         const double* z,     // numObs x 1
                         const double* x,     // numObs * numPredictors
                         size_t numObservations,
                         size_t numPredictors,
                         const bool* z_test, // numTestObs x 1
                         const double* x_test, // numTestObs x (numPredictors + 1)
                         size_t numTestObservations,
                         const double* zetaYs, // numZetaY x 1
                         const double* zetaZs, // numZetaZ x 1
                         size_t numZetaY,
                         size_t numZetaZ,
                         double theta,        // prior prob of U = 1
                         EstimandType estimand,
                         ConfounderModelType confounderModel,
                         size_t numSimsPerCell,
                         size_t numInitialBurnIn,
                         size_t numCellSwitchBurnIn,
                         size_t numTreeSamplesToThin,
                         size_t numThreads,
                         double* estimates,      // numZetaY x numZetaZ x numSimsPerCell
                         double* standardErrors, // numZetaY x numZetaZ
                         bool verbose)
  {
    Control control = { estimand, confounderModel, numSimsPerCell, numInitialBurnIn, numCellSwitchBurnIn,
                        numTreeSamplesToThin, numThreads, theta, verbose, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    lookupBARTFunctions(control);
    
    Data data(y, z, x, numObservations, numPredictors,
              z_test, x_test, numTestObservations);
    
    
    size_t numCells = numZetaY * numZetaZ;
    GridCell* gridCells = new GridCell[numCells];
    
    size_t cellNumber = 0;
    for (size_t i = 0; i < numZetaY; ++i) {
      // snake our way through
      if (i % 2 == 0) {
        for (size_t j = 0; j < numZetaZ; ++j) {
          gridCells[cellNumber].zetaY = zetaYs[i];
          gridCells[cellNumber].zetaZ = zetaZs[j];
          gridCells[cellNumber].offset = i + numZetaY * j;
          gridCells[cellNumber].cellNumber = cellNumber;
          ++cellNumber;
        }
      } else {
        // will wrap around
        for (size_t j = numZetaZ - 1; j < numZetaZ; --j) {
          gridCells[cellNumber].zetaY = zetaYs[i];
          gridCells[cellNumber].zetaZ = zetaZs[j];
          gridCells[cellNumber].offset = i + numZetaY * j;
          gridCells[cellNumber].cellNumber = cellNumber;
          ++cellNumber;
        }
      }
    }
    
#ifdef HAVE_GETTIMEOFDAY
    struct timeval startTime;
    struct timeval endTime;
#else
    time_t startTime;
    time_t endTime;
#endif
    
    
    if (control.numThreads <= 1) {
      Scratch scratch(control, data);
      
      ThreadData threadData = { gridCells, numCells, numCells, &control, &data, &scratch, estimates, standardErrors };
      
#ifdef HAVE_GETTIMEOFDAY
      gettimeofday(&startTime, NULL);
#else
      startTime = time(NULL);
#endif
      
      sensitivityAnalysisTask(&threadData);
    
#ifdef HAVE_GETTIMEOFDAY
      gettimeofday(&endTime, NULL);
#else
      endTime = time(NULL);
#endif
    } else {
      ext_mt_manager_t threadManager;
      ext_mt_create(&threadManager, control.numThreads);
      
      size_t numCellsPerThread;
      size_t offByOneIndex;
      
      ext_mt_getNumThreadsForJob(threadManager, numCells, 0, NULL, &numCellsPerThread, &offByOneIndex);
      
      
      ThreadData* threadData = ext_stackAllocate(control.numThreads, ThreadData);
      void** threadDataPtrs  = ext_stackAllocate(control.numThreads, void*);
      for (size_t i = 0; i < offByOneIndex; ++i) {
        threadData[i].gridCells         = gridCells + i * numCellsPerThread;
        threadData[i].numGridCells      = numCellsPerThread;
        threadData[i].totalNumGridCells = numCells;
        threadData[i].control = &control;
        threadData[i].data    = &data;
        threadData[i].scratch = new Scratch(control, data);
        threadData[i].estimates = estimates;
        threadData[i].standardErrors = standardErrors;
        threadDataPtrs[i] = threadData + i;
      }
      for (size_t i = offByOneIndex; i < numThreads; ++i) {
        threadData[i].gridCells         = gridCells + offByOneIndex * numCellsPerThread + (i - offByOneIndex) * (numCellsPerThread - 1);
        threadData[i].numGridCells      = numCellsPerThread - 1;
        threadData[i].totalNumGridCells = numCells;
        threadData[i].control = &control;
        threadData[i].data    = &data;
        threadData[i].scratch = new Scratch(control, data);
        threadData[i].estimates = estimates;
        threadData[i].standardErrors = standardErrors;
        threadDataPtrs[i] = threadData + i;
      }
      
#ifdef HAVE_GETTIMEOFDAY
      gettimeofday(&startTime, NULL);
#else
      startTime = time(NULL);
#endif
      
      ext_mt_runTasks(threadManager, sensitivityAnalysisTask, threadDataPtrs, control.numThreads);
      
#ifdef HAVE_GETTIMEOFDAY
      gettimeofday(&endTime, NULL);
#else
      endTime = time(NULL);
#endif
      
      for (size_t i = 0; i < control.numThreads; ++i) delete threadData[i].scratch;
      
      ext_stackFree(threadDataPtrs);
      ext_stackFree(threadData);
      
      ext_mt_destroy(threadManager);
    }
    
    if (control.verbose) ext_printf("running time (seconds): %f\n", subtractTimes(endTime, startTime));
    
    delete [] gridCells;
  }
}

namespace {
  using namespace cibart;
  
  void bartCallbackFunction(void* v_callbackData, dbarts::BARTFit& fit, bool,
                            const double* trainingSamples, const double*, double sigma);
  
  struct CallbackData {
    const Control& control;
    const Data& data;
    Scratch& scratch;
    double zetaY, zetaZ;
  };
  
  void sensitivityAnalysisTask(void* v_data)
  {
    ThreadData& threadData(*static_cast<ThreadData*>(v_data));
    
    double* estimates      = threadData.estimates;
    double* standardErrors = threadData.standardErrors;
    const GridCell* gridCells = threadData.gridCells;
    
    Control& control(*threadData.control);
    Data& data(*threadData.data);
    Scratch& scratch(*threadData.scratch);
    
    
    
    sampleConfounders(data, scratch);
    subtractConfounderFromResponse(data, scratch, gridCells[0].zetaY);
    double sigmaEstimate = estimateSigma(data, scratch);
    
    uint32_t* maxNumCuts = new uint32_t[data.numBartPredictors];
    for (size_t i = 0; i < data.numBartPredictors; ++i) maxNumCuts[i] = DEFAULT_BART_MAX_NUM_CUTS;
    
    dbarts::VariableType* variableTypes = new dbarts::VariableType[data.numBartPredictors];
    for (size_t i = 0; i < data.numBartPredictors; ++i) variableTypes[i] = dbarts::ORDINAL;
    
    
    dbarts::Data bartData(scratch.yMinusZetaU, data.bart_x_train, data.x_test, NULL, NULL,
                        data.numObservations, data.numBartPredictors, data.numTestObservations,
                        sigmaEstimate, variableTypes, maxNumCuts);
    
    dbarts::Control bartControl; // use defaults
    bartControl.verbose = false; // otherwise, waaaay too much
    bartControl.numSamples = control.numSimsPerCell;
    bartControl.numBurnIn = control.numInitialBurnIn;
    bartControl.numThreads = 1;
    bartControl.treeThinningRate = control.numTreeSamplesToThin;
    
    dbarts::Model bartModel;
    
    dbarts::CGMPrior* treePrior = ext_stackAllocate(1, dbarts::CGMPrior);
    control.initializeCGMPrior(treePrior, DBARTS_DEFAULT_TREE_PRIOR_BASE, DBARTS_DEFAULT_TREE_PRIOR_POWER);
    
    dbarts::NormalPrior* muPrior = ext_stackAllocate(1, dbarts::NormalPrior);
    control.initializeNormalPrior(muPrior, &bartControl, DBARTS_DEFAULT_NORMAL_PRIOR_K);
    
    dbarts::ChiSquaredPrior* sigmaSqPrior = ext_stackAllocate(1, dbarts::ChiSquaredPrior);
    control.initializeChiSquaredPrior(sigmaSqPrior, DBARTS_DEFAULT_CHISQ_PRIOR_DF, DBARTS_DEFAULT_CHISQ_PRIOR_QUANTILE);
    
    bartModel.treePrior = treePrior;
    bartModel.muPrior = muPrior;
    bartModel.sigmaSqPrior = sigmaSqPrior;
    
    CallbackData callbackData = { control, data, scratch, gridCells[0].zetaY, gridCells[0].zetaZ };
    bartControl.callback = &bartCallbackFunction;
    bartControl.callbackData = (void*) &callbackData;
    
    dbarts::BARTFit* fit = ext_stackAllocate(1, dbarts::BARTFit);
    control.initializeFit(fit, &bartControl, &bartModel, &bartData);
    
    size_t estimateOffset = gridCells[0].offset * control.numSimsPerCell;
    
    dbarts::Results* bartResults = control.runSampler(fit);
    estimateTreatmentEffect(control, data, bartResults->trainingSamples, bartResults->testSamples, estimates + estimateOffset);
    delete bartResults;
    if (control.verbose) ext_printf("Completed %u cell of %u cells.\n", gridCells[0].cellNumber + 1, threadData.totalNumGridCells);
    
    for (size_t i = 1; i < threadData.numGridCells; ++i) {
      // simply change the input BART uses, and make sure our call back function knows
      // what's up
      callbackData.zetaY = gridCells[i].zetaY;
      callbackData.zetaZ = gridCells[i].zetaZ;
      
      subtractConfounderFromResponse(data, scratch, threadData.gridCells[i].zetaY);
      control.setResponse(fit, scratch.yMinusZetaU);
      
      bartResults = control.runSamplerForIterations(fit, control.numCellSwitchBurnIn, control.numSimsPerCell);
      
      estimateTreatmentEffect(control, data, bartResults->trainingSamples, bartResults->testSamples, estimates + estimateOffset);
      standardErrors[gridCells[i].offset] = std::sqrt(ext_computeVariance(estimates + estimateOffset, control.numSimsPerCell, NULL));
      
      delete bartResults;
      if (control.verbose) ext_printf("Completed %u cell of %u cells.\n", gridCells[i].cellNumber + 1, threadData.totalNumGridCells);
    }
    
    control.invalidateFit(fit);
    ext_stackFree(fit);
    
    control.invalidateChiSquaredPrior(sigmaSqPrior);
    ext_stackFree(sigmaSqPrior);
    control.invalidateNormalPrior(muPrior);
    ext_stackFree(muPrior);
    control.invalidateCGMPrior(treePrior);
    ext_stackFree(treePrior);
    
    delete [] maxNumCuts;
    delete [] variableTypes;
  }
  
  // this gets called after BART has new samples for us
  void bartCallbackFunction(void* v_callbackData, dbarts::BARTFit& fit, bool,
                            const double* trainingSamples, const double*, double sigma)
  {
    CallbackData* callbackData = (CallbackData*) v_callbackData;
    
    const Control& control(callbackData->control);
    const Data& data(callbackData->data);
    Scratch& scratch(callbackData->scratch);
    
    updateConfounderProbabilities(control, data, scratch, callbackData->zetaY, callbackData->zetaZ,
                                  trainingSamples, sigma);
    sampleConfounders(data, scratch);
    subtractConfounderFromResponse(data, scratch, callbackData->zetaY);

    control.setResponse(&fit, scratch.yMinusZetaU);
  }
  
  void sampleConfounders(const Data& data, Scratch& scratch)
  {
    for (size_t i = 0; i < data.numObservations; ++i) scratch.u[i] = ext_simulateBernoulli(scratch.p[i]);
  }
  
  void subtractConfounderFromResponse(const Data& data, Scratch& scratch, double zetaY)
  {
    ext_addVectors((const double*) scratch.u, data.numObservations, zetaY, data.y, scratch.yMinusZetaU);
  }
  
  double estimateSigma(const Data& data, Scratch& scratch)
  {
    // perform a standard linear regression
    // we can use the whole X matrix that we allocated consisting of [ 1 X Z ]
    size_t numPredictors = data.numGLMPredictors + 1;
    const double* const& lm_x(data.glm_x);
    
    double* lsSolution = ext_stackAllocate(numPredictors, double);
    double* residuals = ext_stackAllocate(data.numObservations, double);
    char* lsMessage;
    
    int32_t lsResult = ext_findLeastSquaresFit(scratch.yMinusZetaU, data.numObservations, lm_x, numPredictors,
                                               lsSolution, 1.0e-7, residuals, &lsMessage);
    if (lsResult <= 0) ext_throwError("error estimating sigma: %s", lsMessage);
    
    double sumOfSquaredResiduals = ext_sumSquaresOfVectorElements(residuals, data.numObservations);
    
    ext_stackFree(residuals);
    ext_stackFree(lsSolution);
    
    return std::sqrt(sumOfSquaredResiduals / (double) (data.numObservations - data.numBartPredictors));
  }
  
  void estimateTreatmentEffect(const Control& control, const Data& data,
                               const double* trainingSamples, const double* testSamples,
                               double* estimates)
  {
    switch (control.estimand) {
      case ATE:
      {
        double diff;
        for (size_t i = 0; i < control.numSimsPerCell; ++i) {
          double ate = 0.0;
          for (size_t j = 0; j < data.numObservations; ++j) {
            diff = trainingSamples[j] - testSamples[j];
            
            ate += (data.z[j] == 1.0 ? diff : -diff);
          }
          estimates[i] = ate / (double) data.numObservations;
          
          trainingSamples += data.numObservations;
          testSamples     += data.numObservations;
        }
      }
        break;
      case ATT:
      {
        for (size_t i = 0; i < control.numSimsPerCell; ++i) {
          double att = 0.0;
          size_t testIndex = 0;
          for (size_t j = 0; j < data.numObservations; ++j) {
            if (data.z[j] == 0.0) continue;
            
            att += trainingSamples[j] - testSamples[testIndex++];
          }
          estimates[i] = att / (double) data.numTestObservations;
        }
      }
        break;
      case ATC:
      {
        for (size_t i = 0; i < control.numSimsPerCell; ++i) {
          double atc = 0.0;
          size_t testIndex = 0;
          for (size_t j = 0; j < data.numObservations; ++j) {
            if (data.z[j] == 0.0) continue;
            
            atc += testSamples[testIndex++] - trainingSamples[j];
          }
          estimates[i] = atc / (double) data.numTestObservations;
        }
      }
        break;
        
      default:
        ext_throwError("unrecognized estimand type.");
        break;
    }
  }
  
  void updateConfounderProbabilities(const Control& control, const Data& data, Scratch& scratch,
                                     double zetaY, double zetaZ,
                                     const double* yMinusZetaUHat, double sigma)
  {
    double*& glmOffset(scratch.temp_numObs_1);
    ext_scalarMultiplyVector((const double*) scratch.u, data.numObservations, zetaZ, glmOffset);
    
    double* coefficients = ext_stackAllocate(data.numGLMPredictors, double);
    glm_fitGeneralizedLinearModel(data.z, NULL /* sets n = 1 for all obs */, data.numObservations,
                                  data.glm_x, data.numGLMPredictors,
                                  NULL /* weights = 1 */, glmOffset, coefficients,
                                  GLM_FAMILY_BINOMIAL, GLM_LINK_PROBIT, 30,
                                  scratch.glmScratch);
    
    double*& linearPredictor(scratch.temp_numObs_1);
    double zHatForU0, zHatForU1;
    double residForU0, residForU1;
    double pyzu, pyz;
    
    ext_leftMultiplyMatrixAndVector(data.glm_x, data.numObservations, data.numGLMPredictors, coefficients, linearPredictor);
    
    for (size_t i = 0; i < data.numObservations; ++i) {
      zHatForU0 = ext_cumulativeProbabilityOfNormal(linearPredictor[i],         0.0, 1.0);
      zHatForU1 = ext_cumulativeProbabilityOfNormal(linearPredictor[i] + zetaZ, 0.0, 1.0);
      
      residForU0 = data.y[i] - yMinusZetaUHat[i];
      residForU1 = residForU0 - zetaY;
      
      pyzu = ext_densityOfNormal(residForU1, 0.0, sigma) * (data.z[i] == 1.0 ? zHatForU1 : 1.0 - zHatForU1) * control.theta;
      pyz  = ext_densityOfNormal(residForU0, 0.0, sigma) * (data.z[i] == 1.0 ? zHatForU0 : 1.0 - zHatForU0) * (1.0 - control.theta) +
             ext_densityOfNormal(residForU1, 0.0, sigma) * (data.z[i] == 1.0 ? zHatForU1 : 1.0 - zHatForU1) * control.theta;

      scratch.p[i] = pyzu / pyz;
    }
    
    ext_stackFree(coefficients);
  }
  
  Data::Data(const double* y, const double* z, const double* x,
             size_t numObservations, size_t numPredictors, const bool* z_test, const double* x_test,
             size_t numTestObservations) :
    y(y), z(z), x(x), numObservations(numObservations), numPredictors(numPredictors),
    z_test(z_test), x_test(x_test), numTestObservations(numTestObservations),
    bart_x_train(NULL), numBartPredictors(numPredictors + 1),
    glm_x(NULL), numGLMPredictors(numPredictors + 1)
  {
    // we need two different X matrices, built from the one that we are given:
    // [ X Z ] as used in BART and 
    // [ 1 X ] as used in glm
    //
    // since in C matrices are stored by their columns, we can instead allocate
    // [ 1 X Z ]
    // and offset it appropriately
    
    double* temp = new double[numObservations * (numPredictors + 2)];
    
    ext_setVectorToConstant(temp, numObservations, 1.0);
    std::memcpy(temp + numObservations, x, numObservations * numPredictors * sizeof(double));
    std::memcpy(temp + numObservations * (numPredictors + 1), z, numObservations * sizeof(double));
    
    glm_x = temp;
    bart_x_train = temp + numObservations;
  }
  
  Data::~Data() {
    delete [] glm_x;
    glm_x = NULL;
    bart_x_train = NULL;
  }
  
  Scratch::Scratch(const Control& control, const Data& data) : yMinusZetaU(NULL), p(NULL), u(NULL),
    glmScratch(NULL), temp_numObs_1(NULL)
  {
    yMinusZetaU = new double[data.numObservations];
    u = new double[data.numObservations];
    p = new double[data.numObservations];
    ext_setVectorToConstant(p, data.numObservations, control.theta);
    
    glmScratch = new double[glm_getDoubleScratchSize(data.numObservations, data.numGLMPredictors)];
    
    temp_numObs_1 = new double[data.numObservations];
  }
  Scratch::~Scratch()
  {
    delete [] yMinusZetaU; yMinusZetaU = NULL;
    delete [] u; u = NULL;
    delete [] p; p = NULL;
    delete [] glmScratch; glmScratch = NULL;
    delete [] temp_numObs_1; temp_numObs_1 = NULL;
  }
  
  void lookupBARTFunctions(Control& control)
  {
    control.initializeFit             = (void (*)(dbarts::BARTFit*, dbarts::Control*, dbarts::Model*, dbarts::Data*)) R_GetCCallable("dbarts", "initializeFit");
    control.invalidateFit             = (void (*)(dbarts::BARTFit*)) R_GetCCallable("dbarts", "invalidateFit");
    control.runSampler                = (dbarts::Results* (*)(dbarts::BARTFit*)) R_GetCCallable("dbarts", "runSampler");
    control.runSamplerForIterations   = (dbarts::Results* (*)(dbarts::BARTFit*, size_t, size_t)) R_GetCCallable("dbarts", "runSamplerForIterations");
    control.setResponse               = (void (*)(dbarts::BARTFit*, const double*)) R_GetCCallable("dbarts", "setResponse");
    control.initializeCGMPrior        = (void (*)(dbarts::CGMPrior*, double, double)) R_GetCCallable("dbarts", "initializeCGMPriorFromOptions");
    control.invalidateCGMPrior        = (void (*)(dbarts::CGMPrior*)) R_GetCCallable("dbarts", "invalidateCGMPrior");
    control.initializeNormalPrior     = (void (*)(dbarts::NormalPrior*, const dbarts::Control*, double)) R_GetCCallable("dbarts", "initializeNormalPriorFromOptions");
    control.invalidateNormalPrior     = (void (*)(dbarts::NormalPrior*)) R_GetCCallable("dbarts", "invalidateNormalPrior");
    control.initializeChiSquaredPrior = (void (*)(dbarts::ChiSquaredPrior*, double, double)) R_GetCCallable("dbarts", "initializeChiSquaredPriorFromOptions");
    control.invalidateChiSquaredPrior = (void (*)(dbarts::ChiSquaredPrior*)) R_GetCCallable("dbarts", "invalidateChiSquaredPrior");
  }
  
#ifdef HAVE_GETTIMEOFDAY
  double subtractTimes(struct timeval end, struct timeval start) {
    return (1.0e6 * ((double) (end.tv_sec - start.tv_sec)) + (double) (end.tv_usec - start.tv_usec)) / 1.0e6;
  }
#else
  double subtractTimes(time_t end, time_t start) { return (double) (end - start); }
#endif
}
