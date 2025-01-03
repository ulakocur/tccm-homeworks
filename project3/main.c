#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"

int main() {
    	FILE* input_file = fopen("inp.txt", "r");

    	size_t Natoms = read_Natoms(input_file);
    	printf("Number of atoms: %zu\n", Natoms);

    	double** coord = malloc_2d(Natoms, 3);    	
	double* mass = (double*)malloc(Natoms * sizeof(double));
    	read_molecule(input_file, Natoms, coord, mass);
	for (size_t i = 0; i < Natoms; i++) {
		printf("Atom %zu: Coordinates = (%.2f, %.2f, %.2f), Mass = %.2f\n", i + 1, coord[i][0], coord[i][1], coord[i][2], mass[i]);
    }

	double** distances = malloc_2d(Natoms, Natoms);
	compute_distances(Natoms, coord, distances);
	printf("Distances:\n");
	for (size_t i = 0; i < Natoms; i++) {
		for (size_t j = 0; j < Natoms; j++) {
			printf("%.4f ", distances[i][j]);
		}
		printf("\n");
	}

	double sigma = 0.3345;
	double epsilon = 0.0661;
	
	double pot_E = V(epsilon, sigma, Natoms, distances);
	printf("Total potential energy: %.6f J/mol\n", pot_E); 

	double** velocity = malloc_2d(Natoms, 3);
	for (size_t i = 0; i < Natoms; i++) {
		velocity[i][0] = 0.0;
		velocity[i][1] = 0.0;
		velocity[i][2] = 0.0;
	}
	
	double kin_E = T(Natoms, velocity, mass);
	printf("Total kinetic energy: %.6f J/mol\n", kin_E);

	double tot_E = E(epsilon, sigma, Natoms, distances, velocity, mass);
	printf("Total energy: %.6f J/mol\n", tot_E);

	double** acceleration = malloc_2d(Natoms, 3); 
	for (size_t i = 0; i < Natoms; i++) {
                acceleration[i][0] = 0.0;
                acceleration[i][1] = 0.0;
                acceleration[i][2] = 0.0;
        }
	compute_acc(Natoms, coord, mass, distances, acceleration, epsilon, sigma);
	for (size_t i = 0; i < Natoms; i++) {
		printf("Atom %zu: Acceleration = (%.6f, %.6f, %.6f)\n",i+1, acceleration[i][0], acceleration[i][1], acceleration[i][2]);
	}
		
	
//FREEING ALLOCATED MEMORY
    	free_2d(coord);
    	free(mass);
	free_2d(distances);
	free_2d(velocity);
	free_2d(acceleration);
	coord = NULL;
	mass = NULL;
	distances = NULL;
	velocity = NULL;
	acceleration = NULL;

    	fclose(input_file);
    	return 0;
}
