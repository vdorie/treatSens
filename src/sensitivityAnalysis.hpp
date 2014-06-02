#ifndef CIBART_SENSITIVITY_ANALYSIS_HPP
#define CIBART_SENSITIVITY_ANALYSIS_HPP

#include <cstddef>

namespace cibart {
  enum EstimandType {
    ATE,
    ATT,
    ATC
  };
  
  enum ConfounderModelType {
    NORMAL,
    BINOMIAL
  };
  
  void
  fitSensitivityAnalysis(const double* y,     // numObs x 1
                         const double* z,     // numObs x 1
                         const double* x,     // numObs * numCoef, or NULL
                         std::size_t numObservations,
                         std::size_t numPredictors,
                         const bool* z_test, // numTestObs x 1
                         const double* x_test, // numTestObx x numTestCoef or NULL
                         std::size_t numTestObservations,
                         const double* zetaY, // numZetaY x 1
                         const double* zetaZ, // numZetaZ x 1
                         std::size_t numZetaY,
                         std::size_t numZetaZ,
                         double theta,        // prior prob of U = 1
                         EstimandType estimand,
                         ConfounderModelType confounderModel,
                         std::size_t numSimsPerCell,
                         std::size_t numInitialBurnIn,
                         std::size_t numCellSwitchBurnIn,
                         std::size_t numTreeSamplesToThin,
                         std::size_t numThreads,
                         double* estimates,      // numSimsPerCell x numZetaY x numZetaZ
                         double* standardErrors, // numZetaY x numZetaZ
                         bool verbose);
} // namespace cibart

#endif // CIBART_SENSITIVITY_ANALYSIS_HPP
