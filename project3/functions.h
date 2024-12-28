#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double** malloc_2d(size_t m, size_t n);
void free_2d(double** a);
size_t read_Natoms(FILE* input_file);
void read_molecule(FILE* input_file, size_t Natoms, double** coord, double* mass);
void compute_distances(size_t Natoms, double** coord, double** distances);

#endif
