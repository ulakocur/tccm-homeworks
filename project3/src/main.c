#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"


int main(int argc, char* argv[]) {
	// Ensure an input file is provided
	if (argc != 2) {
		printf("Error: Input file needed as the argument (usage: md_simulation [path_to_input_file)]\n");
		exit(-1);
	}
	
	//Open the input file
	FILE* input_file = fopen(argv[1], "r");

	//Read number of atoms and allocate memory for coordinates and masses:
    	size_t Natoms = read_Natoms(input_file);
    	double** coord = malloc_2d(Natoms, 3);    	
    	double* mass = (double*)malloc(Natoms * sizeof(double));
    
	//Read coordinates and masses	
	read_molecule(input_file, Natoms, coord, mass);
	
	//Checking if the mass corresponds to an argon atom
	for(size_t i = 0; i < Natoms; i++) {
		if (mass[i] != 39.948) {
			printf("Error: Properties of atoms other than argon (mass = 39.948) aren't implemented\n");
			exit(-1);
		}
	}
	
	//Close the input file
    	fclose(input_file);

	//Allocating the memory for distances between atoms and computing them
    	double** distances = malloc_2d(Natoms, Natoms);
    	compute_distances(Natoms, coord, distances);

	//Calculate potential energy
    	double pot_E = V(Natoms, distances);

	//Initialize velocity for all atoms to zero
    	double** velocity = malloc_2d(Natoms, 3);
    	for (size_t i = 0; i < Natoms; i++) {
        	for (size_t d = 0; d < 3; d++) {
            		velocity[i][d] = 0.0;
        	}
   	 }

	//Calculate kinetic and total energy
    	double kin_E = T(Natoms, velocity, mass);
    	double tot_E = E(Natoms, distances, velocity, mass);
	
	//Allocating the memory for acceleration and calculating it
    	double** acceleration = malloc_2d(Natoms, 3); 
    	compute_acc(Natoms, coord, mass, distances, acceleration);

	//Set up simulation parameters, tot_steps = total number of steps, dt = time step, M = output frequency 
    	double dt = 0.2; 
    	size_t tot_steps = 1000;
    	size_t M = 1;

	// Open the output files (trajectory.xyz and full.out) and check if they opened correctly
    	FILE* output_file = fopen("trajectory.xyz", "w");
    	if (output_file == NULL) {
        	printf("Error opening output file.\n");
        	exit(-1);
    	}
    	FILE* full_output = fopen("full.out", "w");
    	if (full_output == NULL) {
        	printf("Error opening full output file.\n");
        	exit(-1);
    	}

	//Write initial conditions and energies to full.out
    	fprintf(full_output, "Initial Setup:\n");
	fprintf(full_output, "Energies: Potential = %.6f, Kinetic = %.6f, Total = %.6f\n", pot_E, kin_E, tot_E);
    	fprintf(full_output, "Number of atoms: %zu\n", Natoms);
    	fprintf(full_output, "Atom type: %s\n", ATOM_TYPE);
    	fprintf(full_output, "Sigma: %.4f, Epsilon: %.4f\n", SIGMA, EPSILON);
    	fprintf(full_output, "Coordinates and masses:\n");
    	for (size_t i = 0; i < Natoms; i++) {
        	fprintf(full_output, "Atom %zu (%s): %.5f %.5f %.5f, Mass: %.5f\n", i + 1, ATOM_TYPE, coord[i][0], coord[i][1], coord[i][2], mass[i]);
    	}
	fprintf(full_output, "Distances:\n");
		for (size_t i = 0; i < Natoms; i++) {
                	for (size_t j = i+1; j < Natoms; j++) {
                        	fprintf(full_output, "Atom %zu (%s) - Atoms %zu (%s): %.5f\n",i+1, ATOM_TYPE, j+1, ATOM_TYPE, distances[i][j]);
                        }
                 }
	fprintf(full_output, "\n");

	//Write intial configuration to trajectory.xyz
    	fprintf(output_file, "%zu\n", Natoms);
    	fprintf(output_file, "#Potential energy = %.6f, Kinetic energy = %.6f, Total energy = %.6f \n", pot_E, kin_E, tot_E);	
    	for (size_t i = 0; i < Natoms; i++) {
        	fprintf(output_file, "%s %.5f, %.5f, %.5f\n", ATOM_TYPE,  coord[i][0], coord[i][1], coord[i][2]);
    	}
    
	//Store previous energy for difference calculation
   	 double prev_E = tot_E;

//MAIN SIMULATION LOOP (VERLET ALGORITHM)
	//Update postitions, distances, velocities and accelerations	
    	for (size_t step = 1; step <= tot_steps; step++){
        	update_position(Natoms, coord, velocity, acceleration, distances, dt);
        	compute_distances(Natoms, coord, distances);
        	update_velocity(Natoms, coord, velocity, acceleration, mass, distances, dt);
        	compute_acc(Natoms, coord, mass, distances, acceleration);
        	update_velocity(Natoms, coord, velocity, acceleration, mass, distances, dt);

		//Every M steps update the energies and write the number of atoms, energies and coordinates to trajectory.xyz and atoms, energies, coordinates, velociites and accelerations to full.out
        	if (step % M == 0) {
            		double kin_E = T(Natoms, velocity, mass);
            		double pot_E = V(Natoms, distances);
            		double tot_E = E(Natoms, distances, velocity, mass);
           			double dE = tot_E - prev_E;

            		fprintf(output_file, "%zu\n", Natoms);
            		fprintf(output_file, "#Step %zu: Potential energy = %.6f, Kinetic energy = %.6f, Total energy = %.6f, Energy difference = %.6f \n",step, pot_E, kin_E, tot_E, dE);
            		for (size_t i = 0; i < Natoms; i++) {
                		fprintf(output_file, "%s %.5f, %.5f, %.5f\n", ATOM_TYPE,  coord[i][0], coord[i][1], coord[i][2]);
            		}
            
            		fprintf(full_output, "Step %zu:\n", step);
           		fprintf(full_output, "Energies: Potential = %.6f, Kinetic = %.6f, Total = %.6f, Energy Difference = %.6f\n", pot_E, kin_E, tot_E, dE);
            		fprintf(full_output, "Coordinates:\n");
            		for (size_t i = 0; i < Natoms; i++) {
                		fprintf(full_output, "Atom %zu (%s): %.5f %.5f %.5f\n", i + 1, ATOM_TYPE, coord[i][0], coord[i][1], coord[i][2]);
            		}
			fprintf(full_output, "Distances:\n");
                        for (size_t i = 0; i < Natoms; i++) {
                                for (size_t j = i+1; j < Natoms; j++) {
                                        fprintf(full_output, "Atom %zu (%s) - Atoms %zu (%s): %.5f\n",i+1, ATOM_TYPE, j+1, ATOM_TYPE, distances[i][j]);
                                }
                        }
            		fprintf(full_output, "Velocities:\n");
            		for (size_t i = 0; i < Natoms; i++) {
                		fprintf(full_output, "Atom %zu (%s): %.5f %.5f %.5f\n", i + 1, ATOM_TYPE, velocity[i][0], velocity[i][1], velocity[i][2]);
            		}
            		fprintf(full_output, "Accelerations:\n");
            		for (size_t i = 0; i < Natoms; i++) {
                		fprintf(full_output, "Atom %zu (%s): %.5f %.5f %.5f\n", i + 1, ATOM_TYPE, acceleration[i][0], acceleration[i][1], acceleration[i][2]);
            		}
            		fprintf(full_output, "\n");
	            	
			//Update previous energy for next step
			prev_E = tot_E;
        	}
    	}

//Close the output files
    	fclose(output_file);
    	fclose(full_output);


// Free the allocated memory
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

//Simulation complete message
	printf("MD simulation complete.\n");
    	return 0;
}
