#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"

int main() {
    	FILE* input_file = fopen("inp.txt", "r");

    	size_t Natoms = read_Natoms(input_file);

    	double** coord = malloc_2d(Natoms, 3);    	
	double* mass = (double*)malloc(Natoms * sizeof(double));
    	read_molecule(input_file, Natoms, coord, mass);

	fclose(input_file);

	double** distances = malloc_2d(Natoms, Natoms);
	compute_distances(Natoms, coord, distances);

	double pot_E = V(Natoms, distances);

	double** velocity = malloc_2d(Natoms, 3);
	for (size_t i = 0; i < Natoms; i++) {
		for (size_t d = 0; d < 3; d++) {
			velocity[i][d] = 0.0;
		}
	}

	
	double kin_E = T(Natoms, velocity, mass);

	double tot_E = E(Natoms, distances, velocity, mass);

	double** acceleration = malloc_2d(Natoms, 3); 
	for (size_t i = 0; i < Natoms; i++) {
		for (size_t d = 0; d < 3; d++) {
                	acceleration[i][d] = 0.0;
        	}
	}

	compute_acc(Natoms, coord, mass, distances, acceleration);

	double dt = 0.2;
        size_t tot_steps = 1000;
        size_t M = 10;

	FILE* output_file = fopen("trajectory.xyz", "w");
	if (output_file == NULL) {
		printf("Error opening output file .\n");
		exit(-1);
	}
	
	fprintf(output_file, "%zu\n", Natoms);
	fprintf(output_file, "#Potential energy = %.6f, Kinetic energy = %.6f, Total energy = %.6f \n", pot_E, kin_E, tot_E);	
	for (size_t i = 0; i < Natoms; i++) {
        	fprintf(output_file, "%s %.5f, %.5f, %.5f\n", ATOM_TYPE,  coord[i][0], coord[i][1], coord[i][2]);
	}
	
	double prev_E = tot_E;

	for (size_t step = 1; step < tot_steps; step++){
		update_position(Natoms, coord, velocity, acceleration, distances, dt);
		compute_distances(Natoms, coord, distances);
		update_velocity(Natoms, coord, velocity, acceleration, mass, distances, dt);
		compute_acc(Natoms, coord, mass, distances, acceleration);
		update_velocity(Natoms, coord, velocity, acceleration, mass, distances, dt);
			
		if (step % M == 0) {
			double kin_E = T(Natoms, velocity, mass);
			double pot_E = V(Natoms, distances);
			double tot_E = E(Natoms, distances, velocity, mass);
			double dE = tot_E - prev_E;
		
			fprintf(output_file, "%zu\n", Natoms);
			fprintf(output_file, "#Potential energy = %.6f, Kinetic energy = %.6f, Total energy = %.6f, Energy difference = %.6f \n", pot_E, kin_E, tot_E, dE);
			for (size_t i = 0; i < Natoms; i++) {
                		fprintf(output_file, "%s %.5f, %.5f, %.5f\n", ATOM_TYPE,  coord[i][0], coord[i][1], coord[i][2]);
	        		}
			prev_E = tot_E;
		}
	}	
	fclose(output_file);
	
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

    	return 0;
}
