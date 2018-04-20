#ifndef RNG_H
#define RNG_H

#include "pcg/pcg_random.hpp"

extern pcg32 g_rng;

// Seed RNG using user-provided value.
//   pcg32::state_type is uint64_t.
void setRngSeed(pcg32::state_type seed);

// Return random number in range [0, 1)
double randomNumber();

#endif