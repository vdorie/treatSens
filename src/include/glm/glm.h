#ifndef GLM_GLM_H
#define GLM_GLM_H

#ifdef __cplusplus
#  include <cstddef>
#  define glm_size_t std::size_t
#else
#  include <stddef.h>
#  define glm_size_t size_t
#endif

#include <stdint.h>  // uint32_t

// y_i is a double vector of responses, possibly successes out of n_i trials
// numObs is the total number of i
// x is numObs x numPreds
//
// IT IS ASSUMED THAT ALL INPUTS ARE KOSHER
// if stuff segfaults, probably on you
// specifically:
//   matrices exist and are of appropriate dimension
//   all values of y are between 0 and 1
//   all prior weights are non-zero
//
// if n == NULL, assumed to be 1
// if weights == NULL, also assumed to be 1
// if scratch == NULL, is new'd and deleted

// feel free to implement others of these
typedef enum {
  GLM_FAMILY_BINOMIAL
} glm_family_t;

typedef enum {
  GLM_LINK_PROBIT
} glm_link_t;

#ifdef __cplusplus
extern "C" {
#endif

  void glm_fitGeneralizedLinearModel(const double* y, const double* n, glm_size_t numObs,
                                     const double* x, glm_size_t numCoefs,
                                     const double* weights, const double* offset,
                                     double* beta,
                                     glm_family_t family, glm_link_t link, uint32_t maxIterations,
                                     double* scratch);
  
  glm_size_t glm_getDoubleScratchSize(glm_size_t numObs, glm_size_t numCoefs);
  
#ifdef __cplusplus
}
#endif

#endif // GLM_GLM_H
