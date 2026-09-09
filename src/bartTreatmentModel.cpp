#include "config.hpp"

#include "bartTreatmentModel.hpp"

#include <cstddef> // size_t

#include <external/random.h>
#include <external/stats.h>

#define DBARTS_USE_STUBS
#include <dbarts/dbarts.h>

using std::size_t;

namespace {
  using cibart::TreatmentModel;
  using cibart::BARTTreatmentModel;

  void* createScratch(TreatmentModel* restrict model, ext_rng* restrict generator, const double* restrict x, size_t numObservations, size_t numPredictors, const double* restrict z);
  void destroyScratch(TreatmentModel* model, void* scratch);
  void updateParameters(TreatmentModel* restrict model, void* restrict scratch, const double* restrict offset);
  void getConditionalProbabilities(TreatmentModel* restrict model, void* restrict scratch, double zetaZ, double* restrict probU0, double* restrict probU1);

  // per-analysis state for the propensity sampler; the sampler continues its
  // chain across updateParameters calls (one sweep each)
  struct Scratch {
    dbarts_sampler* fit; // borrowed from the model
    const double* z; // treatment response, borrowed
    double* zHat;    // latest tree-only latent fit
    size_t numObservations;
  };
}

namespace cibart {
  BARTTreatmentModel::BARTTreatmentModel(dbarts_sampler* fit) :
    fit(fit)
  {
    predictorsIncludeIntercept = false;
    includesLatentVariables = false;
    this->createScratch = &::createScratch;
    this->destroyScratch = &::destroyScratch;
    this->updateParameters = &::updateParameters;
    this->getConditionalProbabilities = &::getConditionalProbabilities;
    this->updateLatentVariables = NULL;
  }

  BARTTreatmentModel::~BARTTreatmentModel() { }
}

namespace {
  void* createScratch(TreatmentModel* restrict modelPtr, ext_rng* restrict, const double* restrict, size_t numObservations, size_t, const double* restrict z)
  {
    BARTTreatmentModel& model(*static_cast<BARTTreatmentModel*>(modelPtr));

    Scratch* scratch = new Scratch;
    scratch->numObservations = numObservations;
    scratch->z = z;
    scratch->zHat = new double[numObservations];

    // the probit sampler for the binary propensity model, created R-side (the
    // engine seeded its own chain RNG from R's stream there)
    scratch->fit = model.fit;
    dbarts_sampler_setNumThreads(scratch->fit, 1);
    dbarts_sampler_setVerbose(scratch->fit, 0, 100);

    return scratch;
  }

  void destroyScratch(TreatmentModel*, void* scratchPtr)
  {
    Scratch* scratch = static_cast<Scratch*>(scratchPtr);
    if (scratch != NULL) {
      // releases the engine early; the R object that owns the handle outlives
      // this call and is dropped by its caller right after
      if (scratch->fit != NULL) dbarts_sampler_destroy(scratch->fit);
      delete [] scratch->zHat;
      delete scratch;
    }
  }

  void updateParameters(TreatmentModel* restrict, void* restrict scratchPtr, const double* restrict offset)
  {
    Scratch* scratch = static_cast<Scratch*>(scratchPtr);
    size_t numObservations = scratch->numObservations;

    dbarts_sampler_setOffset(scratch->fit, offset, 1);

    dbarts_results results = {};
    results.structSize = sizeof(results);
    results.train = scratch->zHat;
    dbarts_sampler_run(scratch->fit, 0, 1, &results);

    // a fit made with an offset carries the offset added in; subtract it to
    // recover the tree-only latent prediction
    for (size_t i = 0; i < numObservations; ++i) scratch->zHat[i] -= offset[i];
  }

  void getConditionalProbabilities(TreatmentModel* restrict, void* restrict scratchPtr, double zetaZ, double* restrict probZForU0, double* restrict probZForU1)
  {
    Scratch* scratch = static_cast<Scratch*>(scratchPtr);
    const double* restrict z = scratch->z;
    const double* restrict zHat = scratch->zHat;
    size_t numObservations = scratch->numObservations;

    for (size_t i = 0; i < numObservations; ++i) {
      double zHatU0 = ext_cumulativeProbabilityOfNormal(zHat[i], 0.0, 1.0);
      double zHatU1 = ext_cumulativeProbabilityOfNormal(zHat[i] + zetaZ, 0.0, 1.0);

      probZForU0[i] = (z[i] == 1.0 ? zHatU0 : 1.0 - zHatU0);
      probZForU1[i] = (z[i] == 1.0 ? zHatU1 : 1.0 - zHatU1);
    }
  }
}
