#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ATOM_TYPE "Ar"
#define SIGMA 0.3345
#define EPSILON 0.0661

double** malloc_2d(size_t m, size_t n);
void free_2d(double** a);
size_t read_Natoms(FILE* input_file);
void read_molecule(FILE* input_file, size_t Natoms, double** coord, double* mass);
void compute_distances(size_t Natoms, double** coord, double** distances);
double V(size_t Natoms, double** distance);
double T(size_t Natoms, double** velocity, double* mass);
double E(size_t Natoms, double** distance, double** velocity, double* mass);
void compute_acc(size_t Natoms, double** coord, double* mass, double** distance, double** acceleration);
void update_position(size_t Natoms, double** coord, double** velocity, double** acceleration,double** distances, double dt);
void update_velocity(size_t Natoms, double** coord, double** velocity, double** acceleration, double* mass, double** distances, double dt);

#endif


