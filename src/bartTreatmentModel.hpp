#ifndef CIBART_BART_TREATMENT_MODEL_HPP
#define CIBART_BART_TREATMENT_MODEL_HPP

#include "treatmentModel.hpp"

#include <cstddef> // size_t

// this unusual set of double declarations solves a rather obscure warning on Solaris
extern "C" {
  typedef void* (*C_voidPtrFunction)(void);
  typedef C_voidPtrFunction (*C_voidPtrFunctionLookup)(const char* _namespace, const char* name);
}

namespace cibart {
  struct BARTTreatmentModelFunctionTable;
  
  typedef C_voidPtrFunctionLookup voidPtrFunctionLookup;
  
  struct BARTTreatmentModel : TreatmentModel {
    std::size_t numTrees;
    std::size_t numThin;
    
    double nodePriorParameter;
    
    BARTTreatmentModelFunctionTable* functionTable;
    
    BARTTreatmentModel(voidPtrFunctionLookup lookup, std::size_t numTrees, std::size_t numThin, double nodePriorParameter);
    ~BARTTreatmentModel();
  };
}

#endif // CIBART_BART_TREATMENT_MODEL_HPP
