#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"

int main() {
    	FILE* input_file = fopen("inp.txt", "r");

    	size_t num_atoms = read_Natoms(input_file);
    	printf("Number of atoms: %zu\n", num_atoms);

    	double** coord = malloc_2d(num_atoms, 3);    	
	double* mass = (double*)malloc(num_atoms * sizeof(double));
    	read_molecule(input_file, num_atoms, coord, mass);
	for (size_t i = 0; i < num_atoms; i++) {
		printf("Atom %zu: Coordinates = (%.2f, %.2f, %.2f), Mass = %.2f\n", i + 1, coord[i][0], coord[i][1], coord[i][2], mass[i]);
    }

	double** distances = malloc_2d(num_atoms, num_atoms);
	compute_distances(num_atoms, coord, distances);
	printf("Distances:\n");
	for (size_t i = 0; i < num_atoms; i++) {
		for (size_t j = 0; j < num_atoms; j++) {
			printf("%.2f ", distances[i][j]);
		}
		printf("\n");
	}

	double sigma = 0.3345;
	double epsilon = 0.0661;
	
	double pot_E = V(epsilon, sigma, num_atoms, distances);
	printf("Total potential energy: %.6f J/mol\n", pot_E); 



//FREEING ALLOCATED MEMORY
    	free_2d(coord);
    	free(mass);
	free_2d(distances);
	coord = NULL;
	mass = NULL;
	distances = NULL;

    	fclose(input_file);
    	return 0;
}
