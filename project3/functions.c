#include "functions.h"

//NUMBER OF ATOMS
size_t read_Natoms(FILE* input_file) {
    if (input_file == NULL) {
        printf("Error opening input file\n");
        exit(-1);
    }

    size_t number_of_atoms;
    int read_number = fscanf(input_file, "%zu", &number_of_atoms);
    if (read_number != 1) {
        printf("Error: Not a valid number of atoms\n");
        exit(-1);
    }

    return number_of_atoms;
}

//READING COORDINATES AND MASS
void read_molecule(FILE* input_file, size_t Natoms, double** coord, double* mass) {
    for (size_t i = 0; i < Natoms; i++) {
        int read_number = fscanf(input_file, "%lf %lf %lf %lf", &coord[i][0], &coord[i][1], &coord[i][2], &mass[i]);
        if (read_number != 4) {
            printf("Error: Couldn't read mass and coordinates\n");
            exit(-1);
        }
    }
}
