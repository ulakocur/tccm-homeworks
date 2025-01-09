#include "functions.h"

//ALLOCATE 2D ARRAYS
// Function to allocate memory for a 2D array, where m = number of rows and n = numbers of columns
double** malloc_2d(size_t m, size_t n) {
	double** a = malloc(m*sizeof(double*));
	if (a == NULL) {
		return NULL;
	}
	a[0] = malloc(n*m*sizeof(double));
	if (a[0] == NULL) {
		free(a);
		return NULL;
	}
	for (size_t i=1 ; i<m ; i++) {
		a[i] = a[i-1]+n;
	}
	return a;
}

//FREE 2D ARRAYS
//Free previously allocated memeory for a 2D array
void free_2d(double** a) {
	free(a[0]);
	a[0] = NULL;
	free(a);
}

//READING NUMBER OF ATOMS
//Check if the input file was opened correctly, exit the program with an error code otherwise
size_t read_Natoms(FILE* input_file) {
    if (input_file == NULL) {
        printf("Error opening input file\n");
        exit(-1);
    }

//Read the number of atoms from the first line
    size_t number_of_atoms;
    int read_number = fscanf(input_file, "%zu", &number_of_atoms);
    if (read_number != 1) {
        printf("Error: Not a valid number of atoms\n");
        exit(-1);
    }

    return number_of_atoms;
}

//READING COORDINATES AND MASS
//Read coordinates and mass from the input file
void read_molecule(FILE* input_file, size_t Natoms, double** coord, double* mass) {
    for (size_t i = 0; i < Natoms; i++) {
        int read_number = fscanf(input_file, "%lf %lf %lf %lf", &coord[i][0], &coord[i][1], &coord[i][2], &mass[i]);
//Exit the program if the number of columns in the input file is different than 4
	if (read_number != 4) {
            printf("Error: Couldn't read mass and coordinates\n");
            exit(-1);
        }
    }
}

//COMPUTE DISTANCES BETWEEN ATOMS
//Calculates the distances between all pairs of atoms
void compute_distances(size_t Natoms, double** coord, double** distances) {
	for (size_t i = 0; i < Natoms; i++) {
		for (size_t j = 0; j < Natoms; j++) {
			if (i == j) {
				distances[i][j] = 0.0;
			}
			else {
				double dx = coord[i][0] - coord[j][0];
				double dy = coord[i][1] - coord[j][1];
				double dz = coord[i][2] - coord[j][2];
				distances[i][j] = sqrt(dx * dx + dy * dy + dz * dz);
			}
		}
	}
}

//COMPUTE POTENTIAL ENERGY
//Computes the Lennard-Jones potential energy for all pairs of atoms, checks if all atoms have distinct positions
double V(size_t Natoms, double** distance) {
	double pot_E = 0.0;
	
	for (size_t i = 0; i < Natoms; i++) {
		for (size_t j = i +1; j < Natoms; j++) {
			double r = distance[i][j];
			if (r == 0.0) {
				printf("Error: Multiple atoms in the same position\n");
        			exit(-1);
			}
			double s_r = SIGMA /r;
			double s_r_6 = pow(s_r, 6);
			double s_r_12 = pow(s_r, 12);
	
			double V_lj = 4 * EPSILON * (s_r_12 - s_r_6);

			pot_E += V_lj;
		}
	}
	return pot_E;
}

//COMPUTE KINETIC ENERGY
//Calculates the kinetic energy of each atom
double T(size_t Natoms, double** velocity, double* mass) {
	double kin_E = 0.0;

	for (size_t i = 0; i < Natoms; i++) {
		double v_sq = pow(velocity[i][0], 2) + pow(velocity[i][1], 2) + pow(velocity[i][2], 2);

		kin_E += 0.5 * mass[i] * v_sq;
	}
	return kin_E;
}

//COMPUTE TOTAL ENERGY
//Adds kinetic and potential energy
double E(size_t Natoms, double** distance, double** velocity, double* mass) {
	double pot_E = V(Natoms, distance);
	double kin_E = T(Natoms, velocity, mass);
	double tot_E = kin_E + pot_E;
	return tot_E;
}

//COMPUTE ACCELERATION
//Calculates acceleration vectors for each atom
void compute_acc(size_t Natoms, double** coord, double* mass, double** distance, double** acceleration) {
        for (size_t i = 0; i < Natoms; i++) {
		for (size_t d =0; d < 3; d++) {
			acceleration[i][d] = 0.0;
		}

                for (size_t j = 0; j < Natoms; j++) {
			if (i == j) continue; //No self-interaction

                        double r = distance[i][j];

                        double s_r = SIGMA /r;
                        double s_r_6 = pow(s_r, 6);
                        double s_r_12 = pow(s_r, 12);

                        double U = 24 * EPSILON * (s_r_6 - 2 * s_r_12)/r;

			for (size_t d = 0; d < 3; d++) {
                        	acceleration[i][d] -= U * (coord[i][d] - coord[j][d]) / (mass[i] * r);
			}
		}
        }
}

//UPDATING THE POSITIONS FOR THE VERLET ALGORITHM
void update_position(size_t Natoms, double** coord, double** velocity, double** acceleration, double** distances, double dt) {
	for (size_t i = 0; i < Natoms; i++) {
		for (size_t d = 0; d < 3; d++) {
			coord[i][d] += velocity[i][d] * dt + 0.5 * acceleration[i][d] * pow(dt, 2);
		}
	}
}

//UPDATING THE VELOCITIES FOR THE VERLET ALGORITHM
void update_velocity(size_t Natoms, double** coord, double** velocity, double** acceleration, double* mass, double** distances, double dt) {
	for (size_t i = 0; i < Natoms; i++) {
		for (size_t d = 0; d < 3; d++) {
			velocity[i][d] += 0.5 * acceleration[i][d] * dt;
		}
	}
}
