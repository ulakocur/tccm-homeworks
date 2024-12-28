#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <stdio.h>
#include <stdlib.h>

size_t read_Natoms(FILE* input_file);
void read_molecule(FILE* input_file, size_t Natoms, double** coord, double* mass);

#endif
