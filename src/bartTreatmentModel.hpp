#ifndef CIBART_BART_TREATMENT_MODEL_HPP
#define CIBART_BART_TREATMENT_MODEL_HPP

#include "treatmentModel.hpp"

#include <cstddef> // size_t

#include <external/Rinternals.h> // SEXP

namespace cibart {
  // Optional propensity model: a probit BART fit of Z on X, driven through the
  // flat C API (dbarts.h). The fully-resolved dbarts spec triple is built in R
  // and borrowed here (kept alive R-side for the analysis's duration).
  struct BARTTreatmentModel : TreatmentModel {
    SEXP controlExpr;
    SEXP modelExpr;
    SEXP dataExpr;

    BARTTreatmentModel(SEXP controlExpr, SEXP modelExpr, SEXP dataExpr);
    ~BARTTreatmentModel();
  };
}

#endif // CIBART_BART_TREATMENT_MODEL_HPP
