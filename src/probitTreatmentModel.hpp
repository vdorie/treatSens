#ifndef CIBART_PROBIT_TREATMENT_MODEL_HPP
#define CIBART_PROBIT_TREATMENT_MODEL_HPP

#include "treatmentModel.hpp"

#include <cstddef> // size_t

namespace cibart {
  enum ProbitPriorType {
    PROBIT_PRIOR_FLAT,
    PROBIT_PRIOR_NORMAL,
    PROBIT_PRIOR_STUDENT_T
  };
  
  // keep these POD or else destructor mayhem
  struct ProbitPrior {
  };
  
  struct ProbitNormalPrior : ProbitPrior {
    const double* scale;
  };
  
  struct ProbitStudentTPrior : ProbitPrior {
    const double* scale;
    double dof; // name because *someone* decided defining 'df' as a macro was a good idea
  };
  
  struct ProbitTreatmentModel : TreatmentModel {
    ProbitPriorType priorType;
    const ProbitPrior* prior;
    
    // if FLAT, prior should just be NULL
    ProbitTreatmentModel(ProbitPriorType priorType, const ProbitPrior* prior);
  };
}

#endif // CIBART_PROBIT_TREATMENT_MODEL_HPP
