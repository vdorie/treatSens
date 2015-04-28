#ifndef CIBART_TREATMENT_MODEL_HPP
#define CIBART_TREATMENT_MODEL_HPP

/* Since a model isn't a "first class object", and really just some code we want to inject
   at a particular point in an algorithm, we don't need full-blown virtual functions and 
   an inheritence structure. A treatment model then is just a struct with a few sensible
   flags at the front, some functions in places we can find, and whatever else desired
   after that. Consider them C classes and not C++. */

#include <external/stddef.h> // size_t, restrict

struct ext_rng;

namespace cibart {
  struct TreatmentModel {
    bool predictorsIncludeIntercept;
    bool includesLatentVariables;
    
    void* (*createScratch)(TreatmentModel* restrict model, ext_rng* restrict generator, const double* restrict x, std::size_t numObservations, std::size_t numPredictors,
                           const double* restrict z);
    void (*destroyScratch)(TreatmentModel* model, void* scratch);
    // offset is a temp in cibart's memory; make a copy if persistence is needed
    void (*updateParameters)(TreatmentModel* restrict model, void* restrict scratch, const double* restrict offset);
    void (*getConditionalProbabilities)(TreatmentModel* restrict model, void* restrict scratch, double zetaZ, double* restrict probU0, double* restrict probU1);
    void (*updateLatentVariables)(TreatmentModel* restrict model, void* restrict scratch, const double* restrict offset);
  };
}

#endif // CIBART_TREATMENT_MODEL_HPP
