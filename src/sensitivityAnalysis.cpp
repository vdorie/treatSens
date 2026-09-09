// necessary to get UINT_LEAST32_MAX
#define __STDC_LIMIT_MACROS 1

#include "config.hpp"
#include "sensitivityAnalysis.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring> // memcpy

#if !defined(HAVE_SYS_TIME_H) && defined(HAVE_GETTIMEOFDAY)
#  undef HAVE_GETTIMEOFDAY
#endif
#ifdef HAVE_SYS_TIME_H
#  include <sys/time.h> // gettimeofday
#else
#  include <time.h>
#endif

#include <misc/alloca.h>
#include <misc/linearAlgebra.h>
#include <misc/stats.h>

#include <external/io.h>
#include <external/linearAlgebra.h>
#include <external/random.h>
#include <external/stats.h>

#include "treatmentModel.hpp"

#define DBARTS_USE_STUBS
#include <dbarts/dbarts.h>

#if __cplusplus < 201112L
#  if defined(_WIN64) || SIZEOF_SIZE_T == 8
#    define SIZE_T_SPECIFIER "%lu"
#  else
#    define SIZE_T_SPECIFIER "%u"
#  endif
#else
#  define SIZE_T_SPECIFIER "%zu"
#endif

using std::size_t;
using std::uint32_t;

namespace {
  using namespace cibart;

  // The plain configuration; the classic engine's in-C++ Control/Model/Data and
  // the function-pointer table are gone with the old ABI - the outcome sampler
  // is created R-side and driven here through the flat C API (dbarts.h).
  struct Control {
    EstimandType estimand;
    TreatmentModel& treatmentModel;

    size_t numSimsPerCell;
    size_t numInitialBurnIn;
    size_t numCellSwitchBurnIn;
    size_t numTreeSamplesToThin;

    double theta; // prior probability of U == 1

    bool verbose;
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
    const double* x_train; // root matrix [ 1 X Z ], used for the sigma estimate

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

    double* temp_numObs_1;
    double* temp_numObs_2;

    // landing buffers for discarded (burn-in) sweeps
    double* discardTrain;
    double* discardTest;

    ext_rng* rng;

    Scratch(const Control& control, const Data& data, uint_least32_t rngSeed);
    ~Scratch();
  };

  struct GridCell {
    double zetaY;
    double zetaZ;
    size_t offset;
    size_t cellNumber;
  };

  // forward declarations
  void sampleConfounders(const Data& data, Scratch& scratch);
  void subtractConfounderFromResponse(const Data& data, Scratch& scratch, double zetaY);
  double estimateSigma(const Data& data, Scratch& scratch);
  void updateTreatmentModelParameters(const Control& control, const Data& data, Scratch& scratch, double zetaZ);
  void updateTreatmentModelLatentVariables(const Control& control, const Data& data, Scratch& scratch, double zetaZ);
  void updateConfounderProbabilities(const Control& control, const Data& data, Scratch& scratch, double zetaY, double zetaZ,
                                     const double* yMinusZetaUHat, double sigma);
  void estimateTreatmentEffect(const Control& control, const Data& data,
                               const double* trainingSamples, const double* testSamples, double* estimates);

  // one Gibbs sweep of the confounder model in the driver: update the treatment
  // model, recompute the confounder probabilities from the current BART fit,
  // draw new confounders and rebuild the shifted response. Mirrors the old
  // post-sweep BART callback, less the setResponse the run loop now issues.
  void conditionalGibbsUpdate(const Control& control, const Data& data, Scratch& scratch,
                              double zetaY, double zetaZ, const double* trainingSamples, double sigma)
  {
    updateTreatmentModelParameters(control, data, scratch, zetaZ);
    updateConfounderProbabilities(control, data, scratch, zetaY, zetaZ, trainingSamples, sigma);
    sampleConfounders(data, scratch);
    subtractConfounderFromResponse(data, scratch, zetaY);
    if (control.treatmentModel.includesLatentVariables)
      updateTreatmentModelLatentVariables(control, data, scratch, zetaZ);
  }

#ifdef HAVE_GETTIMEOFDAY
  double subtractTimes(struct timeval end, struct timeval start);
#else
  double subtractTimes(time_t end, time_t start);
#endif

  // Runs one grid cell: numBurnIn discarded sweeps then numSimsPerCell recorded
  // sweeps, each an explicit dbarts_sampler_run(fit, 0, 1, .) followed by the
  // driver's Gibbs update and a setResponse for the next sweep.
  void runGridCell(const Control& control, const Data& data, Scratch& scratch,
                   dbarts_sampler* fit, const GridCell& cell, size_t numBurnIn,
                   double* trainStore, double* testStore,
                   double* estimates, double* standardErrors)
  {
    const double zetaY = cell.zetaY;
    const double zetaZ = cell.zetaZ;
    const size_t numObservations = data.numObservations;
    const size_t numTestObservations = data.numTestObservations;

    dbarts_results results = {};
    results.structSize = sizeof(results);
    double sigmaDraw;

    for (size_t b = 0; b < numBurnIn; ++b) {
      results.sigma = &sigmaDraw;
      results.train = scratch.discardTrain;
      results.test  = numTestObservations > 0 ? scratch.discardTest : NULL;

      dbarts_sampler_run(fit, 0, 1, &results);

      conditionalGibbsUpdate(control, data, scratch, zetaY, zetaZ, results.train, sigmaDraw);
      dbarts_sampler_setResponse(fit, scratch.yMinusZetaU, 1);
    }

    for (size_t s = 0; s < control.numSimsPerCell; ++s) {
      results.sigma = &sigmaDraw;
      results.train = trainStore + s * numObservations;
      results.test  = numTestObservations > 0 ? testStore + s * numTestObservations : NULL;

      dbarts_sampler_run(fit, 0, 1, &results);

      conditionalGibbsUpdate(control, data, scratch, zetaY, zetaZ, results.train, sigmaDraw);
      dbarts_sampler_setResponse(fit, scratch.yMinusZetaU, 1);
    }

    size_t estimateOffset = cell.offset * control.numSimsPerCell;
    estimateTreatmentEffect(control, data, trainStore, testStore, estimates + estimateOffset);
    standardErrors[cell.offset] = std::sqrt(misc_computeVariance(estimates + estimateOffset, control.numSimsPerCell, NULL));
  }
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
                         dbarts_sampler* fit,
                         uint_least32_t rngSeed,
                         double* estimates,      // numZetaY x numZetaZ x numSimsPerCell
                         double* standardErrors, // numZetaY x numZetaZ
                         bool verbose)
  {
    Control control = { estimand, treatmentModel, numSimsPerCell, numInitialBurnIn,
                        numCellSwitchBurnIn, numTreeSamplesToThin, theta, verbose };

    Data data(y, z, x, numObservations, numPredictors, x_test, numTestObservations);

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

    // the confounder Gibbs draws run on a dedicated, R-seeded generator (the
    // classic setRNGState path is gone); the BART engines draw from their own
    // chain RNGs, seeded from R's stream at creation
    Scratch scratch(control, data, rngSeed);

    // the outcome sampler: gaussian family, created R-side once and continued
    // across cells via setResponse; force single-threaded, inline execution
    dbarts_sampler_setNumThreads(fit, 1);
    dbarts_sampler_setVerbose(fit, 0, 100);

    double* trainStore = new double[numObservations * numSimsPerCell];
    double* testStore = numTestObservations > 0 ? new double[numTestObservations * numSimsPerCell] : NULL;

#ifdef HAVE_GETTIMEOFDAY
    struct timeval startTime, endTime;
    gettimeofday(&startTime, NULL);
#else
    time_t startTime = time(NULL), endTime;
#endif

    // cell 0: seed the confounders, response and sigma before the first sweep
    sampleConfounders(data, scratch);
    subtractConfounderFromResponse(data, scratch, gridCells[0].zetaY);
    double sigmaEstimate = estimateSigma(data, scratch);
    dbarts_sampler_setSigma(fit, sigmaEstimate);
    dbarts_sampler_setResponse(fit, scratch.yMinusZetaU, 1);

    runGridCell(control, data, scratch, fit, gridCells[0], numInitialBurnIn,
                trainStore, testStore, estimates, standardErrors);
    if (verbose) {
      ext_printf("Completed cell " SIZE_T_SPECIFIER " of " SIZE_T_SPECIFIER " cells.\n", gridCells[0].cellNumber + 1, numCells);
      ext_fflush_stdout();
    }

    for (size_t i = 1; i < numCells; ++i) {
      // simply change the response BART sees; the Gibbs update inside the cell
      // conditions on the new zetaY / zetaZ
      subtractConfounderFromResponse(data, scratch, gridCells[i].zetaY);
      dbarts_sampler_setResponse(fit, scratch.yMinusZetaU, 1);

      runGridCell(control, data, scratch, fit, gridCells[i], numCellSwitchBurnIn,
                  trainStore, testStore, estimates, standardErrors);
      if (verbose) {
        ext_printf("Completed cell " SIZE_T_SPECIFIER " of " SIZE_T_SPECIFIER " cells.\n", gridCells[i].cellNumber + 1, numCells);
        ext_fflush_stdout();
      }
    }

#ifdef HAVE_GETTIMEOFDAY
    gettimeofday(&endTime, NULL);
#else
    endTime = time(NULL);
#endif
    if (verbose) ext_printf("running time (seconds): %f\n", subtractTimes(endTime, startTime));

    delete [] testStore;
    delete [] trainStore;
    dbarts_sampler_destroy(fit);
    delete [] gridCells;
  }
}

namespace {
  using namespace cibart;

  void sampleConfounders(const Data& data, Scratch& scratch)
  {
    for (size_t i = 0; i < data.numObservations; ++i) scratch.u[i] = static_cast<double>(ext_rng_simulateBernoulli(scratch.rng, scratch.p[i]));
  }

  void subtractConfounderFromResponse(const Data& data, Scratch& scratch, double zetaY)
  {
    misc_addVectors(static_cast<const double*>(scratch.u), data.numObservations, -zetaY, data.y, scratch.yMinusZetaU);
  }

  double estimateSigma(const Data& data, Scratch& scratch)
  {
    // perform a standard linear regression
    // we can use the whole X matrix that we allocated consisting of [ 1 X Z ]
    size_t numPredictors = data.numPredictors + 2;
    const double* const& lm_x(data.x_train);

    double* lsSolution = misc_stackAllocate(numPredictors, double);
    double* residuals = misc_stackAllocate(data.numObservations, double);
    char* lsMessage;

    int32_t lsResult = ext_findLeastSquaresFit(scratch.yMinusZetaU, data.numObservations, lm_x, numPredictors,
                                               lsSolution, 1.0e-7, residuals, &lsMessage);
    if (lsResult <= 0) ext_throwError("error estimating sigma: %s", lsMessage);

    double sumOfSquaredResiduals = ext_sumSquaresOfVectorElements(residuals, data.numObservations);

    misc_stackFree(residuals);
    misc_stackFree(lsSolution);

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
    misc_scalarMultiplyVector(const_cast<const double*>(scratch.u), data.numObservations, zetaZ, offset);

    control.treatmentModel.updateParameters(&control.treatmentModel, scratch.treatmentScratch, offset);
  }

  void updateTreatmentModelLatentVariables(const Control& control, const Data& data, Scratch& scratch, double zetaZ)
  {
    double*& offset(scratch.temp_numObs_1);
    misc_scalarMultiplyVector(const_cast<const double*>(scratch.u), data.numObservations, zetaZ, offset);

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
    x_test(x_test), numTestObservations(numTestObservations), x_train(NULL)
  {
    // create matrix [ 1 X Z ]; the sigma estimate uses [ 1 X Z ], the outcome
    // BART's predictors [ X Z ] and test matrix are built R-side into its spec
    double* x_temp = new double[numObservations * (numPredictors + 2)];
    misc_setVectorToConstant(x_temp, numObservations, 1.0);
    std::memcpy(x_temp + numObservations, x, numObservations * numPredictors * sizeof(double));
    std::memcpy(x_temp + numObservations * (numPredictors + 1), z, numObservations * sizeof(double));

    x_train = x_temp;
  }

  Data::~Data() {
    delete [] x_train;
    x_train = NULL;
  }

  Scratch::Scratch(const Control& control, const Data& data, uint_least32_t rngSeed) :
    yMinusZetaU(NULL), p(NULL), u(NULL),
    treatmentModel(control.treatmentModel), treatmentScratch(NULL),
    temp_numObs_1(NULL), temp_numObs_2(NULL), discardTrain(NULL), discardTest(NULL)
  {
    // a standalone generator (default algorithm + standard-normal), reseeded
    // deterministically from R's stream so set.seed governs reproducibility
    rng = ext_rng_createDefault(false);
    ext_rng_setSeed(rng, rngSeed);

    yMinusZetaU = new double[data.numObservations];
    u = new double[data.numObservations];
    p = new double[data.numObservations];
    misc_setVectorToConstant(p, data.numObservations, control.theta);

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

    discardTrain = new double[data.numObservations];
    discardTest = data.numTestObservations > 0 ? new double[data.numTestObservations] : NULL;
  }

  Scratch::~Scratch()
  {
    delete [] discardTest; discardTest = NULL;
    delete [] discardTrain; discardTrain = NULL;

    delete [] temp_numObs_2; temp_numObs_2 = NULL;
    delete [] temp_numObs_1; temp_numObs_1 = NULL;

    treatmentModel.destroyScratch(&treatmentModel, treatmentScratch);
    treatmentScratch = NULL;

    delete [] p; p = NULL;
    delete [] u; u = NULL;
    delete [] yMinusZetaU; yMinusZetaU = NULL;

    ext_rng_destroy(rng);
  }

#ifdef HAVE_GETTIMEOFDAY
  double subtractTimes(struct timeval end, struct timeval start) {
    return (1.0e6 * static_cast<double>(end.tv_sec - start.tv_sec) + static_cast<double>(end.tv_usec - start.tv_usec)) / 1.0e6;
  }
#else
  double subtractTimes(time_t end, time_t start) { return static_cast<double>(end - start); }
#endif
}
