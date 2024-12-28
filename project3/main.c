#include <stdio.h>
#include <stdlib.h>
#include "functions.h"

int main() {
    FILE* input_file = fopen("inp.txt", "r");

    size_t num_atoms = read_Natoms(input_file);

    printf("Number of atoms: %zu\n", num_atoms);

    double** coord = (double**)malloc(num_atoms * sizeof(double*));
    for (size_t i = 0; i < num_atoms; i++) {
        coord[i] = (double*)malloc(3 * sizeof(double));
    }
    double* mass = (double*)malloc(num_atoms * sizeof(double));

    read_molecule(input_file, num_atoms, coord, mass);

    for (size_t i = 0; i < num_atoms; i++) {
        printf("Atom %zu: Coordinates = (%.2f, %.2f, %.2f), Mass = %.2f\n", i + 1, coord[i][0], coord[i][1], coord[i][2], mass[i]);
    }

//FREEING ALLOCATED MEMORY
    for (size_t i = 0; i < num_atoms; i++) {
        free(coord[i]);
    }
    free(coord);
    free(mass);

    fclose(input_file);
    return 0;
}
