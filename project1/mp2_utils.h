// mp2_utils.h
#ifndef MP2_UTILS_H
#define MP2_UTILS_H

#include <stdint.h> // Required for uint64_t

// Function declarations
uint64_t encode_indices(int i, int j, int k, int l);
void decode_key(uint64_t key, int* i, int* j, int* k, int* l);

#endif // MP2_UTILS_H

