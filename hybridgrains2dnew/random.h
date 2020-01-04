/*
 * uniformly distributed random number generation with
 * twice faster version of Mersenne Twister.
 *
 * $Id: random.h,v 1.4 2003/06/18 11:35:15 syoyo Exp $
 */

#ifndef RANDOM_H
#define RANDOM_H

//#include "vector.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

extern void seedMT(uint32_t seed);
extern float randomMT(void);
extern uint32_t randomMT_i(void);

// extern void      random_uniform_vector(ri_vector_t *dst);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
