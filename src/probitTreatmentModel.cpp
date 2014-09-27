#include "probitTreatmentModel.hpp"

#include <cstring>

#include <external/stddef.h>
#include <external/stats.h>
#include <external/linearAlgebra.h>

#include <external/io.h>

namespace {
  struct Scratch;
  
  void* createScratch(cibart::TreatmentModel* restrict model, const double* x, size_t numObservations, size_t numPredictors, const double* z);
  void destroyScratch(void* scratch);
  void updateParameters(cibart::TreatmentModel* restrict model, void* restrict scratch, const double* restrict offset);
  void getConditionalProbabilities(cibart::TreatmentModel* restrict model, void* restrict scratch, double zetaZ, double* restrict probU0, double* restrict probU1);
  
  void sampleCoefficients(cibart::ProbitTreatmentModel& restrict model, Scratch& restrict scratch, const double* restrict offset);
  void sampleLatents(cibart::ProbitTreatmentModel& restrict model, Scratch& restrict scratch, const double* restrict offset);
  
  void getPosteriorCovarianceFactor(cibart::ProbitTreatmentModel& model, const double* xCrossproduct, size_t numPredictors, double* factor);
  void updatePosteriorCovarianceFactor(cibart::ProbitTreatmentModel& model, Scratch& scratch);
}

namespace cibart {
  ProbitTreatmentModel::ProbitTreatmentModel(ProbitPriorType priorType, const ProbitPrior* prior) :
    priorType(priorType), prior(prior)
  {
    predictorsIncludeIntercept = true;
    this->createScratch = &::createScratch;
    this->destroyScratch = &::destroyScratch;
    this->updateParameters = &::updateParameters;
    this->getConditionalProbabilities = &::getConditionalProbabilities;
  }
}

namespace {
  struct Scratch {
    const double* x;
    const double* z;
    size_t numObservations;
    size_t numPredictors;
    double* latents;
    double* coefficients;
    double sigma_sq; // used to create t or cauchy prior
    
    // denote the coef as 'beta' and the latent variable as 'v'
    // beta | v ~ N(cov X' (v - offset), cov)
    // cov = (X'X + Sig_beta^{-1})^{-1},
    // Let LL' = X'X + Sig_beta^{-1}, then beta | v ~ N(L^-T L^-1 * X' (v - mu), L^-T L^-1),
    // if b ~ N(0, I), L^-T (u + L^-1 X'(v - mu)) has desired dist
    double* posteriorCovarianceInverseRightFactor; // this is L^T
    double* xCrossproduct; // cached for t dists
  };

  void* createScratch(cibart::TreatmentModel* restrict modelPtr, const double* restrict x, size_t numObservations, size_t numPredictors, const double* restrict z)
  {
    Scratch* scratch = new Scratch;
    cibart::ProbitTreatmentModel& model(*((cibart::ProbitTreatmentModel*) modelPtr));
    
    ext_printf("probit treatment model:\n");
    if (model.priorType == cibart::PROBIT_PRIOR_STUDENT_T) {
      cibart::ProbitStudentTPrior& prior(*((cibart::ProbitStudentTPrior *) model.prior));
      ext_printf("  family = t\n");
      ext_printf("  df = %f\n", prior.dof);
      ext_printf("  scale = %f", prior.scale[0]);
      for (size_t i = 1; i < numPredictors; ++i) ext_printf(", %f", prior.scale[i]);
      ext_printf("\n");
    } else if (model.priorType == cibart::PROBIT_PRIOR_NORMAL) {
      cibart::ProbitNormalPrior& prior(*((cibart::ProbitNormalPrior *) model.prior));
      ext_printf("  family = normal\n");
      ext_printf("  scale = %f", prior.scale[0]);
      for (size_t i = 1; i < numPredictors; ++i) ext_printf(", %f", prior.scale[i]);
      ext_printf("\n");
    } else {
      ext_printf("  family = flat\n");
    }
    R_FlushConsole();
    
    scratch->x = x;
    scratch->numObservations = numObservations;
    scratch->numPredictors   = numPredictors;
    scratch->z = z;
    
    scratch->latents      = new double[numObservations];
    scratch->coefficients = new double[numPredictors];
    scratch->sigma_sq = 1.0;
    
    scratch->posteriorCovarianceInverseRightFactor = new double[numPredictors * numPredictors];
    scratch->xCrossproduct = NULL;
    
    double* xCrossproduct = scratch->posteriorCovarianceInverseRightFactor;
    
    if (model.priorType == cibart::PROBIT_PRIOR_STUDENT_T) {
      // this only needs to be cached for T priors, and is otherwise computed where the factor goes
      scratch->xCrossproduct = new double[numPredictors * numPredictors];
      xCrossproduct = scratch->xCrossproduct;
    }
    
    ext_getSingleMatrixCrossproduct(x, numObservations, numPredictors, xCrossproduct,
                                    false, EXT_TRIANGLE_TYPE_UPPER);
    
    getPosteriorCovarianceFactor(model, xCrossproduct, numPredictors, scratch->posteriorCovarianceInverseRightFactor);
    
    
    // use some default values for latents
    for (size_t i = 0; i < numObservations; ++i) {
      scratch->latents[i] = (z[i] == 1.0 ? 0.5 : -0.5);
    }
    
    return scratch;
  }
  
  void destroyScratch(void* scratchPtr) {
    Scratch* scratch = (Scratch*) scratchPtr;
    
    if (scratch != NULL) {
      delete [] scratch->xCrossproduct;
      delete [] scratch->posteriorCovarianceInverseRightFactor;
      delete [] scratch->coefficients;
      delete [] scratch->latents;
      delete scratch;
    }
  }
  
  void updateParameters(cibart::TreatmentModel* restrict modelPtr, void* restrict scratchPtr, const double* restrict offset)
  {
    cibart::ProbitTreatmentModel& model(*((cibart::ProbitTreatmentModel*) modelPtr));
    Scratch& restrict scratch(*(Scratch*) scratchPtr);
    
    sampleCoefficients(model, scratch, offset);
    sampleLatents(model, scratch, offset);
  }
  
  void getConditionalProbabilities(cibart::TreatmentModel* restrict, void* restrict scratchPtr, double zetaZ, double* restrict probZForU0, double* restrict probZForU1)
  {
    Scratch& restrict scratch(*(Scratch*) scratchPtr);
    
    double* restrict& temp(probZForU0);
    ext_leftMultiplyMatrixAndVector(scratch.x, scratch.numObservations, scratch.numPredictors,
                                    scratch.coefficients, temp);
    
    for (size_t i = 0; i < scratch.numObservations; ++i) {
      double zHatU0 = ext_cumulativeProbabilityOfNormal(temp[i], 0.0, 1.0);
      double zHatU1 = ext_cumulativeProbabilityOfNormal(temp[i] + zetaZ, 0.0, 1.0);
      
      probZForU0[i] = (scratch.z[i] == 1.0 ? zHatU0 : 1.0 - zHatU0);
      probZForU1[i] = (scratch.z[i] == 1.0 ? zHatU1 : 1.0 - zHatU1);
    }
  }
  
  void sampleCoefficients(cibart::ProbitTreatmentModel& restrict model, Scratch& restrict scratch, const double* restrict offset)
  {
    if (model.priorType == cibart::PROBIT_PRIOR_STUDENT_T) {
      cibart::ProbitStudentTPrior& prior(*((cibart::ProbitStudentTPrior*) model.prior));
      
      // sample sigma_sq first:
      double sigmaScale = 1.0;
      for (size_t i = 0; i < scratch.numPredictors; ++i) {
        sigmaScale += (scratch.coefficients[i] / prior.scale[i]) * (scratch.coefficients[i] / prior.scale[i]);
      }
      sigmaScale /= prior.dof;
      
      scratch.sigma_sq = sigmaScale / ext_simulateChiSquared(prior.dof + (double) scratch.numPredictors);
      
      updatePosteriorCovarianceFactor(model, scratch);
    }
    
    // v - zetaZ * u
    ext_addVectorsInPlace(offset, scratch.numObservations, -1.0, scratch.latents);
    // X'(v - zetaZ * u)
    ext_multiplyMatrixIntoVector(scratch.x, scratch.numObservations, scratch.numPredictors, true, scratch.latents, scratch.coefficients);
    
    // L^{-1} X'(v - zetaZ * u)
    ext_solveTriangularSystemInPlace(scratch.posteriorCovarianceInverseRightFactor, scratch.numPredictors, true, EXT_TRIANGLE_TYPE_UPPER,
                                     scratch.coefficients, 1);

    // b + L^{-1} X'(v - zetaZ * u), b ~ N(0, 1)
    for (size_t i = 0; i < scratch.numPredictors; ++i) {
      scratch.coefficients[i] += ext_simulateStandardNormal();
    }
    
    // beta = L^{-T} (b + L^{-1} X'(v - zetaZ * u))
    ext_solveTriangularSystemInPlace(scratch.posteriorCovarianceInverseRightFactor, scratch.numPredictors, false, EXT_TRIANGLE_TYPE_UPPER,
                                     scratch.coefficients, 1);
  }
  
  void sampleLatents(cibart::ProbitTreatmentModel& restrict model, Scratch& restrict scratch, const double* restrict offset)
  {
    // X beta
    ext_multiplyMatrixIntoVector(scratch.x, scratch.numObservations, scratch.numPredictors, false, scratch.coefficients, scratch.latents);
    // X beta + zetaZ * u
    ext_addVectorsInPlace(offset, scratch.numObservations, 1.0, scratch.latents);
    
    // ~ N(X beta + zetaZ * u, 1)^{+/-}
    for (size_t i = 0; i < scratch.numObservations; ++i) {
      scratch.latents[i] = (scratch.z[i] == 0.0 ?
                            ext_simulateLowerTruncatedNormalScale1(scratch.latents[i], 0.0) : 
                            ext_simulateUpperTruncatedNormalScale1(scratch.latents[i], 0.0));
    }
  }
  
  void getPosteriorCovarianceFactor(cibart::ProbitTreatmentModel& model, const double* xCrossproduct, size_t numPredictors, double* factor)
  {    
    switch (model.priorType) {
      case cibart::PROBIT_PRIOR_STUDENT_T:
      {
        cibart::ProbitStudentTPrior& prior(*((cibart::ProbitStudentTPrior*) model.prior));
        
        std::memcpy(factor, xCrossproduct, numPredictors * numPredictors * sizeof(double));
        
        // X'X + Sig_beta^{-1}, but we only do diagonal Sig_beta 
        for (size_t i = 0; i < numPredictors; ++i) {
          factor[i * (1 + numPredictors)] += 1.0 / (prior.scale[i] * prior.scale[i] * prior.dof);
        }
      }
      break;
      case cibart::PROBIT_PRIOR_NORMAL:
      {
        cibart::ProbitNormalPrior& prior(*((cibart::ProbitNormalPrior*) model.prior));

        for (size_t i = 0; i < numPredictors; ++i) {
          factor[i * (1 + numPredictors)] += 1.0 / (prior.scale[i] * prior.scale[i]);
        }
      }
      break;
      default:
      break;
    }
    
    ext_getSymmetricPositiveDefiniteTriangularFactorizationInPlace(factor, numPredictors, EXT_TRIANGLE_TYPE_UPPER);
    
    // fill in lower left as zeros
    for (size_t col = 0; col < numPredictors - 1; ++col) {
      for (size_t row = col + 1; row < numPredictors; ++row) factor[row + col * numPredictors] = 0.0;
    }
  }
  
  // only applicable to T prior
  void updatePosteriorCovarianceFactor(cibart::ProbitTreatmentModel& model, Scratch& scratch)
  {
    const double* restrict xCrossproduct = scratch.xCrossproduct;
    double* restrict factor = scratch.posteriorCovarianceInverseRightFactor;
    size_t numPredictors = scratch.numPredictors;
    
    cibart::ProbitStudentTPrior& prior(*((cibart::ProbitStudentTPrior*) model.prior));
    
    std::memcpy(factor, xCrossproduct, numPredictors * numPredictors * sizeof(double));
    
    // X'X + Sig_beta^{-1}, but we only do diagonal Sig_beta 
    for (size_t i = 0; i < numPredictors; ++i) {
      factor[i * (1 + numPredictors)] += 1.0 / (prior.scale[i] * prior.scale[i] * prior.dof * scratch.sigma_sq);
    }
    
    ext_getSymmetricPositiveDefiniteTriangularFactorizationInPlace(factor, numPredictors, EXT_TRIANGLE_TYPE_UPPER);
    
    // fill in lower left as zeros
    for (size_t col = 0; col < numPredictors - 1; ++col) {
      for (size_t row = col + 1; row < numPredictors; ++row) factor[row + col * numPredictors] = 0.0;
    }
  }
}
