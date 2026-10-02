#include "config.hpp"

#include "probitEMTreatmentModel.hpp"

#include <misc/stddef.h>
#include <external/stats.h>
#include <external/random.h>
#include <external/linearAlgebra.h>
#include <external/Rinternals.h> // R_alloc
#include <glm/glm.h>

#include <external/io.h>

namespace {
  void* createScratch(cibart::TreatmentModel* restrict model, ext_rng* restrict generator, const double* restrict x, size_t numObservations, size_t numPredictors, const double* restrict z);
  void destroyScratch(cibart::TreatmentModel*, void* scratch);
  void updateParameters(cibart::TreatmentModel* restrict model, void* restrict scratch, const double* restrict offset);
  void getConditionalProbabilities(cibart::TreatmentModel* restrict model, void* restrict scratch, double zetaZ, double* restrict probU0, double* restrict probU1);
}

namespace cibart {
  ProbitEMTreatmentModel::ProbitEMTreatmentModel(size_t maxIterations) :
    maxIterations(maxIterations)
  {
    predictorsIncludeIntercept = true;
    includesLatentVariables = false;
    this->createScratch = &::createScratch;
    this->destroyScratch = &::destroyScratch;
    this->updateParameters = &::updateParameters;
    this->getConditionalProbabilities = &::getConditionalProbabilities;
    updateLatentVariables = NULL;
  }
}

namespace {
  struct Scratch {
    const double* x;
    const double* z;
    size_t numObservations;
    size_t numPredictors;
    double* coefficients;
    double* glmScratch;
  };

  // the scratch and its buffers are R transient storage (R_alloc), so a raise
  // during the analysis - the glm fit's own included - strands none of it;
  // destroyScratch has nothing to free
  void* createScratch(cibart::TreatmentModel* restrict, ext_rng* restrict, const double* restrict x, size_t numObservations, size_t numPredictors, const double* restrict z)
  {
    Scratch* scratch = reinterpret_cast<Scratch*>(R_alloc(1, sizeof(Scratch)));
    
    /* ext_printf("probit EM model:\n");
    cibart::ProbitEMTreatmentModel& model(*((cibart::ProbitEMTreatmentModel*) modelPtr));
    ext_printf("  max iter: %z\n", model.maxIterations); */
    
    scratch->x = x;
    scratch->numObservations = numObservations;
    scratch->numPredictors   = numPredictors;
    scratch->z = z;
    
    scratch->coefficients = reinterpret_cast<double*>(R_alloc(numPredictors, sizeof(double)));
    scratch->glmScratch = reinterpret_cast<double*>(R_alloc(glm_getDoubleScratchSize(numObservations, numPredictors), sizeof(double)));
    
    return scratch;
  }
  
  void destroyScratch(cibart::TreatmentModel*, void*) { }
  
  void updateParameters(cibart::TreatmentModel* restrict modelPtr, void* restrict scratchPtr, const double* restrict offset)
  {
    cibart::ProbitEMTreatmentModel* model = static_cast<cibart::ProbitEMTreatmentModel*>(modelPtr);
    Scratch& restrict scratch(*static_cast<Scratch*>(scratchPtr));
    
    glm_fitGeneralizedLinearModel(scratch.z, NULL /* sets n = 1 for all obs */, scratch.numObservations,
                                  scratch.x, scratch.numPredictors,
                                  NULL /* weights = 1 */, offset, scratch.coefficients,
                                  GLM_FAMILY_BINOMIAL, GLM_LINK_PROBIT, static_cast<uint32_t>(model->maxIterations),
                                  scratch.glmScratch);
  }
  
  void getConditionalProbabilities(cibart::TreatmentModel* restrict, void* restrict scratchPtr, double zetaZ, double* restrict probZForU0, double* restrict probZForU1)
  {
    Scratch& restrict scratch(*static_cast<Scratch*>(scratchPtr));
    
    ext_leftMultiplyMatrixAndVector(scratch.x, scratch.numObservations, scratch.numPredictors,
                                    scratch.coefficients, probZForU0);
    
    for (size_t i = 0; i < scratch.numObservations; ++i) {
      double zHatU0 = ext_cumulativeProbabilityOfNormal(probZForU0[i], 0.0, 1.0);
      double zHatU1 = ext_cumulativeProbabilityOfNormal(probZForU0[i] + zetaZ, 0.0, 1.0);
      
      probZForU0[i] = (scratch.z[i] == 1.0 ? zHatU0 : 1.0 - zHatU0);
      probZForU1[i] = (scratch.z[i] == 1.0 ? zHatU1 : 1.0 - zHatU1);
    }
  }
}
