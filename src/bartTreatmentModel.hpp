#ifndef CIBART_BART_TREATMENT_MODEL_HPP
#define CIBART_BART_TREATMENT_MODEL_HPP

#include "treatmentModel.hpp"

#include <cstddef> // size_t

// dbarts.h's opaque sampler handle, forward declared so this header needs
// neither dbarts.h nor an R header
struct dbarts_sampler_t;

namespace cibart {
  // Optional propensity model: a probit BART fit of Z on X, driven through the
  // flat C API (dbarts.h). dbarts.h declares no creation entry, so the sampler
  // is built in R from its spec triple and the handle is borrowed here (the R
  // object is kept alive for the analysis's duration).
  struct BARTTreatmentModel : TreatmentModel {
    ::dbarts_sampler_t* fit;

    BARTTreatmentModel(::dbarts_sampler_t* fit);
    ~BARTTreatmentModel();
  };
}

#endif // CIBART_BART_TREATMENT_MODEL_HPP
