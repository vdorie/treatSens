// necessary to get UINT_LEAST32_MAX
#define __STDC_LIMIT_MACROS 1

#include "config.hpp"
#include "sensitivityAnalysis.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <dbarts/cstdint.hpp>

#if !defined(HAVE_SYS_TIME_H) && defined(HAVE_GETTIMEOFDAY)
#  undef HAVE_GETTIMEOFDAY
#endif
#ifdef HAVE_SYS_TIME_H
#  include <sys/time.h> // gettimeofday
#else
#  include <time.h>
#endif

// #include <R.h>
// #include <Rdefines.h>
// #include <Rinternals.h>
// #include <Rmath.h>
#include <R_ext/Rdynload.h>

// #include <external/Rinternals.h> // SEXP

#include <external/alloca.h>
#include <external/io.h>
#include <external/linearAlgebra.h>
#include <external/random.h>
#include <external/stats.h>
#include <external/thread.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/types.hpp>
#include <dbarts/results.hpp>
#include <dbarts/R_C_interface.hpp>

#include "treatmentModel.hpp"

#define DEFAULT_BART_MAX_NUM_CUTS 100u

using std::size_t;
using std::uint32_t;
using std::uint_least32_t;

namespace {
  using namespace cibart;
  struct Control {
    EstimandType estimand;
    TreatmentModel& treatmentModel;
    
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
    
    const double* x_test;
    size_t numTestObservations;
        
    // various transformations
    const double* x_train; // root matrix
    const double* bart_x_train; // possibly offset
    size_t numBartPredictors;
    
    Data(const double* y, const double* z, const double* x,
         size_t numObservations, size_t numPredictors, const double* x_test,
         size_t numTestObservations);
    ~Data();
  };
  
  struct Scratch {
    double* yMinusZetaU;
    double* p;
    double* u;
    
    TreatmentModel& treatmentModel;
    void* treatmentScratch;
    // void (*destroyTreatmentScratch)(void* treatmentScratch);
    
    double* temp_numObs_1;
    double* temp_numObs_2;
    
    ext_rng* rng;
    
    Scratch(const Control& control, const Data& data);
    ~Scratch();
  };
  
  // forward declarations
  //void getEstimatesForGridPoint(const Control& control, const Data& data, Scratch& scratch, double zetaY, double zetaZ,
  //                              double* estimates);
  void sampleConfounders(const Data& data, Scratch& scratch);
  void subtractConfounderFromResponse(const Data& data, Scratch& scratch, double zetaY);
  double estimateSigma(const Data& data, Scratch& scratch);
  void updateTreatmentModelParameters(const Control& control, const Data& data, Scratch& scratch, double zetaZ);
  void updateTreatmentModelLatentVariables(const Control& control, const Data& data, Scratch& scratch, double zetaZ);
  void updateConfounderProbabilities(const Control& control, const Data& data, Scratch& scratch, double zetaY, double zetaZ,
                                     const double* yMinusZetaUHat, double sigma);
  void estimateTreatmentEffect(const Control& control, const Data& data,
                               const double* trainingSamples, const double* testSamples, double* estimates);
  void lookupBARTFunctions(Control& control);
}

extern "C" {
  static void sensitivityAnalysisTask(void* v_data);
}

namespace {
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
                         const double* x_test, // numTestObs x (numPredictors + 1)
                         size_t numTestObservations,
                         const double* zetaYs, // numZetaY x 1
                         const double* zetaZs, // numZetaZ x 1
                         size_t numZetaY,
                         size_t numZetaZ,
                         double theta,        // prior prob of U = 1
                         EstimandType estimand,
                         TreatmentModel& treatmentModel,
                         size_t numSimsPerCell,
                         size_t numInitialBurnIn,
                         size_t numCellSwitchBurnIn,
                         size_t numTreeSamplesToThin,
                         size_t numThreads,
                         double* estimates,      // numZetaY x numZetaZ x numSimsPerCell
                         double* standardErrors, // numZetaY x numZetaZ
                         bool verbose)
  {
    Control control = { estimand, treatmentModel, numSimsPerCell, numInitialBurnIn, numCellSwitchBurnIn,
                        numTreeSamplesToThin, numThreads, theta, verbose, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
      
    lookupBARTFunctions(control);
    
    Data data(y, z, x, numObservations, numPredictors,
              x_test, numTestObservations);
    
    size_t numCells = numZetaY * numZetaZ;
    if (numCells < numThreads) {
      numThreads = numCells;
      control.numThreads = numCells;
    }
    
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
  void bartCallbackFunctionForLatents(void* v_callbackData, dbarts::BARTFit& fit, bool,
                                      const double* trainingSamples, const double*, double sigma);
  
  struct CallbackData {
    const Control& control;
    const Data& data;
    Scratch& scratch;
    double zetaY, zetaZ;
  };
}

extern "C" {
  static void sensitivityAnalysisTask(void* v_data)
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
    
    
    dbarts::Data bartData(scratch.yMinusZetaU, data.bart_x_train, data.x_test, NULL, NULL, NULL,
                          data.numObservations, data.numBartPredictors, data.numTestObservations,
                          sigmaEstimate, variableTypes, maxNumCuts);
    
    dbarts::Control bartControl; // use defaults
    bartControl.verbose = false; // otherwise, waaaay too much
    bartControl.numSamples = control.numSimsPerCell;
    bartControl.numBurnIn = control.numInitialBurnIn;
    bartControl.numThreads = 1;
    bartControl.treeThinningRate = static_cast<uint32_t>(control.numTreeSamplesToThin);
    bartControl.rng = scratch.rng;
    
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
    bartControl.callback = control.treatmentModel.includesLatentVariables? &bartCallbackFunctionForLatents : &bartCallbackFunction;
    bartControl.callbackData = static_cast<void*>(&callbackData);
    
    dbarts::BARTFit* fit = ext_stackAllocate(1, dbarts::BARTFit);
    control.initializeFit(fit, &bartControl, &bartModel, &bartData);
    
    size_t estimateOffset = gridCells[0].offset * control.numSimsPerCell;
    
    
    dbarts::Results* bartResults = control.runSampler(fit);
    
    estimateTreatmentEffect(control, data, bartResults->trainingSamples, bartResults->testSamples, estimates + estimateOffset);
    standardErrors[gridCells[0].offset] = std::sqrt(ext_computeVariance(estimates + estimateOffset, control.numSimsPerCell, NULL));
    
    delete bartResults;

    if (control.verbose && control.numThreads == 1) {
      ext_printf("Completed cell %u of %u cells.\n", gridCells[0].cellNumber + 1, threadData.totalNumGridCells);
      ext_fflush_stdout();
    }
    
    for (size_t i = 1; i < threadData.numGridCells; ++i) {
      // simply change the input BART uses, and make sure our call back function knows
      // what's up
      callbackData.zetaY = gridCells[i].zetaY;
      callbackData.zetaZ = gridCells[i].zetaZ;
      
      subtractConfounderFromResponse(data, scratch, threadData.gridCells[i].zetaY);
      control.setResponse(fit, scratch.yMinusZetaU);
      
      bartResults = control.runSamplerForIterations(fit, control.numCellSwitchBurnIn, control.numSimsPerCell);
      
      estimateOffset = gridCells[i].offset * control.numSimsPerCell;
      
      estimateTreatmentEffect(control, data, bartResults->trainingSamples, bartResults->testSamples, estimates + estimateOffset);
      standardErrors[gridCells[i].offset] = std::sqrt(ext_computeVariance(estimates + estimateOffset, control.numSimsPerCell, NULL));
      
      delete bartResults;
      if (control.verbose && control.numThreads == 1) {
        ext_printf("Completed cell %u of %u cells.\n", gridCells[i].cellNumber + 1, threadData.totalNumGridCells);
        ext_fflush_stdout();
      }
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
}

namespace {
  // this gets called after BART has new samples for us
  void bartCallbackFunction(void* v_callbackData, dbarts::BARTFit& fit, bool,
                            const double* trainingSamples, const double*, double sigma)
  {
    CallbackData& callbackData(*static_cast<CallbackData*>(v_callbackData));
    
    const Control& control(callbackData.control);
    const Data& data(callbackData.data);
    Scratch& scratch(callbackData.scratch);
        
    updateTreatmentModelParameters(control, data, scratch, callbackData.zetaZ);
    updateConfounderProbabilities(control, data, scratch, callbackData.zetaY, callbackData.zetaZ,
                                  trainingSamples, sigma);
    sampleConfounders(data, scratch);
    subtractConfounderFromResponse(data, scratch, callbackData.zetaY);

    control.setResponse(&fit, scratch.yMinusZetaU);
  }
  
  // this gets called after BART has new samples for us
  void bartCallbackFunctionForLatents(void* v_callbackData, dbarts::BARTFit& fit, bool,
                                      const double* trainingSamples, const double*, double sigma)
  {
    CallbackData& callbackData(*static_cast<CallbackData*>(v_callbackData));
    
    const Control& control(callbackData.control);
    const Data& data(callbackData.data);
    Scratch& scratch(callbackData.scratch);
        
    updateTreatmentModelParameters(control, data, scratch, callbackData.zetaZ);
    updateConfounderProbabilities(control, data, scratch, callbackData.zetaY, callbackData.zetaZ,
                                  trainingSamples, sigma);
    sampleConfounders(data, scratch);
    subtractConfounderFromResponse(data, scratch, callbackData.zetaY);
    updateTreatmentModelLatentVariables(control, data, scratch, callbackData.zetaZ);

    control.setResponse(&fit, scratch.yMinusZetaU);
  }
  
  void sampleConfounders(const Data& data, Scratch& scratch)
  {
    for (size_t i = 0; i < data.numObservations; ++i) scratch.u[i] = static_cast<double>(ext_rng_simulateBernoulli(scratch.rng, scratch.p[i]));
  }
  
  void subtractConfounderFromResponse(const Data& data, Scratch& scratch, double zetaY)
  {
    ext_addVectors(static_cast<const double*>(scratch.u), data.numObservations, -zetaY, data.y, scratch.yMinusZetaU);
  }
  
  double estimateSigma(const Data& data, Scratch& scratch)
  {
    // perform a standard linear regression
    // we can use the whole X matrix that we allocated consisting of [ 1 X Z ]
    size_t numPredictors = data.numPredictors + 2;
    const double* const& lm_x(data.x_train);
    
    double* lsSolution = ext_stackAllocate(numPredictors, double);
    double* residuals = ext_stackAllocate(data.numObservations, double);
    char* lsMessage;
    
    int32_t lsResult = ext_findLeastSquaresFit(scratch.yMinusZetaU, data.numObservations, lm_x, numPredictors,
                                               lsSolution, 1.0e-7, residuals, &lsMessage);
    if (lsResult <= 0) ext_throwError("error estimating sigma: %s", lsMessage);
    
    double sumOfSquaredResiduals = ext_sumSquaresOfVectorElements(residuals, data.numObservations);
            
    ext_stackFree(residuals);
    ext_stackFree(lsSolution);
    
    return std::sqrt(sumOfSquaredResiduals / static_cast<double>(data.numObservations - numPredictors));
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
          estimates[i] = ate / static_cast<double>(data.numObservations);
          
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
          estimates[i] = att / static_cast<double>(data.numTestObservations);
          
          trainingSamples += data.numObservations;
        }
      }
      break;
      case ATC:
      {
        for (size_t i = 0; i < control.numSimsPerCell; ++i) {
          double atc = 0.0;
          size_t testIndex = 0;
          for (size_t j = 0; j < data.numObservations; ++j) {
            if (data.z[j] == 1.0) continue;
            
            atc += testSamples[testIndex++] - trainingSamples[j];
          }
          estimates[i] = atc / static_cast<double>(data.numTestObservations);
          
          trainingSamples += data.numObservations;
        }
      }
      break;
    }
  }
  
  void updateTreatmentModelParameters(const Control& control, const Data& data, Scratch& scratch, double zetaZ)
  {
    double*& offset(scratch.temp_numObs_1);
    ext_scalarMultiplyVector(const_cast<const double*>(scratch.u), data.numObservations, zetaZ, offset);
    
    control.treatmentModel.updateParameters(&control.treatmentModel, scratch.treatmentScratch, offset);
  }
  
  void updateTreatmentModelLatentVariables(const Control& control, const Data& data, Scratch& scratch, double zetaZ)
  {
    double*& offset(scratch.temp_numObs_1);
    ext_scalarMultiplyVector(const_cast<const double*>(scratch.u), data.numObservations, zetaZ, offset);
    
    control.treatmentModel.updateLatentVariables(&control.treatmentModel, scratch.treatmentScratch, offset);
  }

  
  void updateConfounderProbabilities(const Control& control, const Data& data, Scratch& scratch,
                                     double zetaY, double zetaZ,
                                     const double* yMinusZetaUHat, double sigma)
  {
    double*& probZForU0(scratch.temp_numObs_1);
    double*& probZForU1(scratch.temp_numObs_2);
    
    control.treatmentModel.getConditionalProbabilities(&control.treatmentModel, scratch.treatmentScratch, zetaZ, probZForU0, probZForU1);
    
    double probU0, probU1;
    double bartDensityForU0, bartDensityForU1;
    
    for (size_t i = 0; i < data.numObservations; ++i) {
      bartDensityForU0 = ext_densityOfNormal(data.y[i] - yMinusZetaUHat[i], 0.0, sigma);
      bartDensityForU1 = ext_densityOfNormal(data.y[i] - yMinusZetaUHat[i] - zetaY, 0.0, sigma);
      
      probU0 = bartDensityForU0 * probZForU0[i] * (1.0 - control.theta);
      probU1 = bartDensityForU1 * probZForU1[i] * control.theta;

      scratch.p[i] = probU1 / (probU0 + probU1);
    }
  }
  
  Data::Data(const double* y, const double* _z, const double* _x,
             size_t _numObservations, size_t _numPredictors, const double* x_test,
             size_t numTestObservations) :
    y(y), z(_z), x(_x), numObservations(_numObservations), numPredictors(_numPredictors),
    x_test(x_test), numTestObservations(numTestObservations),
    x_train(NULL), bart_x_train(NULL), numBartPredictors(numPredictors + 1)
  {
    // create matrix [ 1 X Z ]; treatmentModel may or may not need [ 1 X ], bart definitely does not
    
    double* x_temp = new double[numObservations * (numPredictors + 2)];
    ext_setVectorToConstant(x_temp, numObservations, 1.0);
    std::memcpy(x_temp + numObservations, x, numObservations * numPredictors * sizeof(double));
    std::memcpy(x_temp + numObservations * (numPredictors + 1), z, numObservations * sizeof(double));
    
    x_train = x_temp;
    bart_x_train = x_temp + numObservations;
  }
  
  Data::~Data() {
    delete [] x_train;
    x_train = NULL;
    bart_x_train = NULL;
  }
  
  Scratch::Scratch(const Control& control, const Data& data) : yMinusZetaU(NULL), p(NULL), u(NULL),
    treatmentModel(control.treatmentModel), treatmentScratch(NULL), temp_numObs_1(NULL), temp_numObs_2(NULL)
  {
    if (control.numThreads > 1) {
      rng = ext_rng_createDefault(false);
    } else {
      rng = ext_rng_createDefault(true);
    }
    
    yMinusZetaU = new double[data.numObservations];
    u = new double[data.numObservations];
    p = new double[data.numObservations];
    ext_setVectorToConstant(p, data.numObservations, control.theta);
    
    // [ 1 X ]
    const double* treatment_x = data.x_train;
    size_t treatmentNumPredictors = data.numPredictors + 1;
    if (!control.treatmentModel.predictorsIncludeIntercept) {
      --treatmentNumPredictors;
      treatment_x += data.numObservations;
    }
    
    treatmentScratch = treatmentModel.createScratch(&treatmentModel, rng, treatment_x, data.numObservations, treatmentNumPredictors, data.z);
    
    temp_numObs_1 = new double[data.numObservations];
    temp_numObs_2 = new double[data.numObservations];
  }
  
  Scratch::~Scratch()
  {
    delete [] temp_numObs_2; temp_numObs_2 = NULL;
    delete [] temp_numObs_1; temp_numObs_1 = NULL;
    
    treatmentModel.destroyScratch(&treatmentModel, treatmentScratch);
    treatmentScratch = NULL;
    
    delete [] p; p = NULL;
    delete [] u; u = NULL;
    delete [] yMinusZetaU; yMinusZetaU = NULL;
    
    ext_rng_destroy(rng);
  }
  
  void lookupBARTFunctions(Control& control)
  {
    control.initializeFit             = reinterpret_cast<void (*)(dbarts::BARTFit*, dbarts::Control*, dbarts::Model*, dbarts::Data*)>(R_GetCCallable("dbarts", "initializeFit"));
    control.invalidateFit             = reinterpret_cast<void (*)(dbarts::BARTFit*)>(R_GetCCallable("dbarts", "invalidateFit"));
    control.runSampler                = reinterpret_cast<dbarts::Results* (*)(dbarts::BARTFit*)>(R_GetCCallable("dbarts", "runSampler"));
    control.runSamplerForIterations   = reinterpret_cast<dbarts::Results* (*)(dbarts::BARTFit*, size_t, size_t)>(R_GetCCallable("dbarts", "runSamplerForIterations"));
    control.setResponse               = reinterpret_cast<void (*)(dbarts::BARTFit*, const double*)>(R_GetCCallable("dbarts", "setResponse"));
    control.initializeCGMPrior        = reinterpret_cast<void (*)(dbarts::CGMPrior*, double, double)>(R_GetCCallable("dbarts", "initializeCGMPriorFromOptions"));
    control.invalidateCGMPrior        = reinterpret_cast<void (*)(dbarts::CGMPrior*)>(R_GetCCallable("dbarts", "invalidateCGMPrior"));
    control.initializeNormalPrior     = reinterpret_cast<void (*)(dbarts::NormalPrior*, const dbarts::Control*, double)>(R_GetCCallable("dbarts", "initializeNormalPriorFromOptions"));
    control.invalidateNormalPrior     = reinterpret_cast<void (*)(dbarts::NormalPrior*)>(R_GetCCallable("dbarts", "invalidateNormalPrior"));
    control.initializeChiSquaredPrior = reinterpret_cast<void (*)(dbarts::ChiSquaredPrior*, double, double)>(R_GetCCallable("dbarts", "initializeChiSquaredPriorFromOptions"));
    control.invalidateChiSquaredPrior = reinterpret_cast<void (*)(dbarts::ChiSquaredPrior*)>(R_GetCCallable("dbarts", "invalidateChiSquaredPrior"));
  }
  
#ifdef HAVE_GETTIMEOFDAY
  double subtractTimes(struct timeval end, struct timeval start) {
    return (1.0e6 * static_cast<double>(end.tv_sec - start.tv_sec) + static_cast<double>(end.tv_usec - start.tv_usec)) / 1.0e6;
  }
#else
  double subtractTimes(time_t end, time_t start) { return static_cast<double>(end - start); }
#endif
}
