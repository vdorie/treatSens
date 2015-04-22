#ifndef CIBART_PROBIT_EM_TREATMENT_MODEL_HPP
#define CIBART_PROBIT_EM_TREATMENT_MODEL_HPP

#include "treatmentModel.hpp"

#include <cstddef> // size_t

namespace cibart {
  struct ProbitEMTreatmentModel : TreatmentModel {
    std::size_t maxIterations;
    
    ProbitEMTreatmentModel(std::size_t maxIterations);
  };
}

#endif // CIBART_PROBIT_EM_TREATMENT_MODEL_HPP
