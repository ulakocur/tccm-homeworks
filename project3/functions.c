#include "functions.h"

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

void free_2d(double** a) {
	free(a[0]);
	a[0] = NULL;
	free(a);
}

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

//COMPUTE DISTANCES BETWEEN ATOMS
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
double V(double epsilon, double sigma, size_t Natoms, double** distance) {
	double pot_E = 0.0;
	
	for (size_t i = 0; i < Natoms; i++) {
		for (size_t j = i +1; j < Natoms; j++) {
			double r = distance[i][j];
			if (r == 0.0) {
				printf("Error: Multiple atoms in the same position\n");
        			exit(-1);
			}
			double s_r = sigma /r;
			double s_r_6 = pow(s_r, 6);
			double s_r_12 = pow(s_r, 12);
	
			double V_lj = 4 * epsilon * (s_r_12 - s_r_6);

			pot_E += V_lj;
		}
	}
	return pot_E;
}

//COMPUTE KINETIC ENERGY
double T(size_t Natoms, double** velocity, double* mass) {
	double kin_E = 0.0;

	for (size_t i = 0; i < Natoms; i++) {
		double v_sq = pow(velocity[i][0], 2) + pow(velocity[i][1], 2) + pow(velocity[i][2], 2);

		kin_E += 0.5 * mass[i] * v_sq;
	}
	return kin_E;
}

//TOTAL ENERGY
double E(double epsilon, double sigma, size_t Natoms, double** distance, double** velocity, double* mass) {
	double pot_E = V(epsilon, sigma, Natoms, distance);
	double kin_E = T(Natoms, velocity, mass);
	double tot_E = kin_E + pot_E;
	return tot_E;
}
