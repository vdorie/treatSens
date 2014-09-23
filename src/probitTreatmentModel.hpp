#ifndef CIBART_PROBIT_TREATMENT_MODEL_HPP
#define CIBART_PROBIT_TREATMENT_MODEL_HPP

#include "treatmentModel.hpp"

#include <cstddef> // size_t

namespace cibart {
  enum ProbitPriorType {
    PROBIT_PRIOR_NORMAL,
    PROBIT_PRIOR_STUDENT_T
  };
  
  struct ProbitPrior {
  };
  
  struct ProbitNormalPrior : ProbitPrior {
    const double* scale;
  };
  
  struct ProbitStudentTPrior : ProbitPrior {
    const double* scale;
    double dof;
  };
  
  struct ProbitTreatmentModel : TreatmentModel {
    ProbitPriorType priorType;
    const ProbitPrior* prior;
    
    ProbitTreatmentModel(ProbitPriorType priorType, const ProbitPrior* prior);
  };
}

#endif // CIBART_PROBIT_TREATMENT_MODEL_HPP
