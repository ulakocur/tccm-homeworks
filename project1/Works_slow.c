#include <stdio.h>
#include <trexio.h>
#include <stdlib.h>
#include <string.h>

// Function to evaluate all terms and permutations
double evaluate_integrals(int* index, double* value, int n_two_elec_int, int n_up);
// Function to evaluate the mp2 corrections
double mp2_energy_correction(int* index, double* value, int n_two_elec_int, int n_up);
double calculate_mp2_energy(int* index, double* value, int n_two_elec_int, int n_up, int n_mo, double* orbital_energies);

int main(int argc, char* argv[]) {
    if (argc != 2) {  // Ensure the user provides 1 argument which is the path to the trexio file
        fprintf(stderr, "Usage: %s <path_to_file>\n", argv[0]);
        return 1;
    }

    const char* file_name = argv[1];  // Get the file path from the command-line argument
    trexio_exit_code rc_open;

    // Open file
    trexio_t* file = trexio_open(file_name, 'r', TREXIO_HDF5, &rc_open);
    if (rc_open != TREXIO_SUCCESS) {
        fprintf(stderr, "Error opening file '%s': %d\n", file_name, rc_open);
        return 1;
    }

    printf("TREXIO file '%s' opened successfully.\n\n", file_name);

    // Read nucleus repulsion energy
    double nuc_rep_energy;
    trexio_exit_code nuc_rep_read = trexio_read_nucleus_repulsion(file, &nuc_rep_energy);
    if (nuc_rep_read != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading nucleus repulsion energy: %d\n", nuc_rep_read);
        trexio_close(file);
        return 1;
    }

    // Read number of molecular orbitals
    int32_t n_orb;
    trexio_exit_code n_orb_read = trexio_read_mo_num(file, &n_orb);
    if (n_orb_read != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading number of orbitals: %d\n", n_orb_read);
        trexio_close(file);
        return 1;
    }

    printf("Number of molecular orbitals: %d\n", n_orb);

    // Read number of non-zero two-electron integrals
    int64_t n_two_elec_int;
    trexio_exit_code n_two_elec_int_read = trexio_read_mo_2e_int_eri_size(file, &n_two_elec_int);
    if (n_two_elec_int_read != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading number of two-electron integrals: %d\n", n_two_elec_int_read);
        trexio_close(file);
        return 1;
    }

    // Read number of spin-up electrons
    int32_t n_up;
    trexio_exit_code n_up_read = trexio_read_electron_up_num(file, &n_up);
    if (n_up_read != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading number of spin-up electrons: %d\n", n_up_read);
        trexio_close(file);
        return 1;
    }
    printf("Number of spin-up electrons: %d\n", n_up);

    // Allocate memory for one-electron integrals
    double* one_e_int_core = malloc(n_orb * n_orb * sizeof(double));
    if (one_e_int_core == NULL) {
        fprintf(stderr, "Memory allocation failed for one-electron integrals.\n");
        trexio_close(file);
        return 1;
    }

    trexio_exit_code one_e_int_core_read = trexio_read_mo_1e_int_core_hamiltonian(file, one_e_int_core);
    if (one_e_int_core_read != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading one-electron integrals.\n");
        free(one_e_int_core);
        trexio_close(file);
        return 1;
    }

    // Allocate memory for two-electron integrals
    int32_t* index = malloc(4 * n_two_elec_int * sizeof(int32_t));
    if (index == NULL) {
        fprintf(stderr, "Memory allocation failed for indices.\n");
        free(one_e_int_core);
        trexio_close(file);
        return 1;
    }

    double* value = malloc(n_two_elec_int * sizeof(double));
    if (value == NULL) {
        fprintf(stderr, "Memory allocation failed for values.\n");
        free(index);
        free(one_e_int_core);
        trexio_close(file);
        return 1;
    }

    trexio_exit_code two_e_int_read = trexio_read_mo_2e_int_eri(file, 0, &n_two_elec_int, index, value);
    if (two_e_int_read != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading two-electron integrals.\n");
        free(index);
        free(value);
        free(one_e_int_core);
        trexio_close(file);
        return 1;
    }

    // Read orbital energies
    double* orbital_energies = malloc(n_orb * sizeof(double));
    if (orbital_energies == NULL) {
        fprintf(stderr, "Memory allocation failed for orbital energies.\n");
        free(index);
        free(value);
        free(one_e_int_core);
        trexio_close(file);
        return 1;
    }

    trexio_exit_code orbital_energies_read = trexio_read_mo_energy(file, orbital_energies);
    if (orbital_energies_read != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading orbital energies.\n");
        free(orbital_energies);
        free(index);
        free(value);
        free(one_e_int_core);
        trexio_close(file);
        return 1;
    }

    // Hartree-Fock energy calculation
    double total_energy = nuc_rep_energy;
    printf("\nStarting energy calculation...\n");
    printf("Nucleus repulsion energy: %f\n", nuc_rep_energy);

    // One-electron energy
    double one_e_energy = 0.0;
    for (int i = 0; i < n_up; i++) {
        one_e_energy += 2.0 * one_e_int_core[i * n_orb + i]; // i * n_orb + i gives only the diagonal which is <i|h|i>
    }
    total_energy += one_e_energy;
    printf("One-Electron Energy: %f\n", one_e_energy);

    // Two-electron energy
    double two_e_energy = evaluate_integrals(index, value, n_two_elec_int, n_up);
    total_energy += two_e_energy;

    printf("Two-Electron Energy: %f\n", two_e_energy);

    // Final total energy
    printf("Total Hartree-Fock Energy: %f\n", total_energy);
    printf("Expected energy: -76.0267987\n\n");

    // MP2 energy correction
    double mp2_energy = calculate_mp2_energy(index, value, n_two_elec_int, n_up, n_orb, orbital_energies);
    printf("MP2 Energy Correction: %f\n", mp2_energy);
    printf("MP2 Correction should be: −0.20395997\n");

    // Add MP2 correction to Hartree-Fock energy
    total_energy += mp2_energy;
    printf("Total Energy with MP2 Correction: %f\n", total_energy);

    // Free memory
    free(one_e_int_core);
    free(index);
    free(value);
    free(orbital_energies);

    // Close the file
    trexio_exit_code rc_close = trexio_close(file);
    if (rc_close != TREXIO_SUCCESS) {
        fprintf(stderr, "Error closing file: %d\n", rc_close);
        return 1;
    }

    printf("TREXIO file '%s' closed successfully.\n", file_name);
    return 0;
}

double evaluate_integrals(int* index, double* value, int n_two_elec_int, int n_up) {
    double two_e_energy = 0.0;

    // Iterate over all two-electron integrals
    for (int m = 0; m < n_two_elec_int; m++) {
        int i = index[4 * m + 0];
        int j = index[4 * m + 1];
        int k = index[4 * m + 2];
        int l = index[4 * m + 3];
        double integral = value[m];

        // Only consider integrals involving occupied orbitals
        if (i < n_up && j < n_up && k < n_up && l < n_up) {
            // Direct term: (ij|ij)
            if (i == k && j == l) { 
                if (i == j) { // this only has a singular permutation which shouldn't be counted double
                    two_e_energy += 2 * integral;
                }
                else {
                    two_e_energy += 2 * 2 * integral; // Direct contribution (is doubled to take the permutations into consideration there are only 2 valid permutations e.g. (1 0 | 1 0) and (0 1 | 0 1) the others reduce to the same)
                }
            }
            // Exchange term: (ij|ji)
            if (i == j && l == k) {
                if (i == k) { // once again only a single permutation which shouldn't be counted double
                    two_e_energy -= integral;
                }
                else {
                    two_e_energy -= 2 * integral; // Direct contribution (Once again there are only 2 valid permutations)
                }
            }
        }
    }

    return two_e_energy;
}

double calculate_mp2_energy(int* index, double* value, int n_two_elec_int, int n_up, int n_mo, double* orbital_energies) {
    double mp2_energy = 0.0;

    // Iterate over all two-electron integrals
    for (int m = 0; m < n_two_elec_int; m++) { // division by 100 to help debug (not print a ton of output)
        int i = index[4 * m + 0];  // Virtual (b, a)
        int j = index[4 * m + 1];
        int k = index[4 * m + 2];  // Occupied (j, i) In reverse since we use permutation (i j|k l) = (l k|j i)
        int l = index[4 * m + 3];
        double integral = value[m];

	// printf("(%d, %d, %d, %d) = %f\n",i ,j ,k ,l ,integral);
        // MP2 only considers (ij|ab) integrals where i, j are occupied, and a, b are virtual, So we need at least 2 variables to be higher or equal than n_up and 2 lower than n_up
        if (((i >= n_up) + (j >= n_up) + (k >= n_up) + (l >= n_up)) == 2) {
	    //printf("Used:(%d, %d, %d, %d) = %f\n",l ,k ,j ,i ,integral);
	    
	    int occupied[2], virtual[2]; // To store the indices
	    int occ_count = 0, virt_count = 0;

	    int indices[4] = {i, j, k, l}; // Group all indices
	    for (int idx = 0; idx < 4; idx++) {
		if (indices[idx] < n_up) {
        	occupied[occ_count++] = indices[idx]; // Add to occupied
    	        } else {
        	virtual[virt_count++] = indices[idx]; // Add to virtual
	        }
	    }
	    //printf("occupied: ");
	    for (int i = 0; i < 2; i++) { // Loop through the array
    		//printf("%d ", occupied[i]);
	    }
	    //printf("virtual: ");
	    for (int i = 0; i < 2; i++) { // Loop through the array
    		//printf("%d ", virtual[i]);
	    }
//printf("\n"); // Add a newline for better readability
	    if ((occupied[0] == i && occupied[1] == k) || (occupied[0] == k && occupied[1] == i) || (occupied[0] == j && occupied[1] == l) || (occupied[0] == l && occupied[1] == j)){
		//printf("This integral should be skipped: (%d, %d ,%d , %d) = %f\n",i, j, k, l, integral);
		continue;
	    }

		double denominator = orbital_energies[occupied[0]] + orbital_energies[occupied[1]] - orbital_energies[virtual[0]] - orbital_energies[virtual[1]];

            // Ensure the denominator is non-zero (Since the degenerate case is much more difficult)
            if (denominator != 0) {
                // Calculate contributions
                double swapped_term = 0.0;

                // Look for the exchange integral <ij|lk> (which is symmetric this can probably be done significantly faster but idk how)
                for (int n = 0; n < n_two_elec_int; n++) {
                    int i2 = index[4 * n + 0];
                    int j2 = index[4 * n + 1];
                    int k2 = index[4 * n + 2];
                    int l2 = index[4 * n + 3];

                    // Check if the integral is <ij|lk> (or a variant from symmetry, but they shouldn't be since they are unique species with either i==j k==l or 4 different entries. Since we search from low numbers to high we should come across it, if I don't find it the integral is actually 0 right?????)
                    if (i2 == i && j2 == j && k2 == l && l2 == k) {
                        swapped_term = value[n];
                        break;
                    }
// I can also just try all permutations? Not necessary I guess but should be foolproof?

		   if (i2 == i && j2 == k && k2 == l && l2 == j) {
                        swapped_term = value[n];
                        break;
                    }
		   if (i2 == l && j2 == k && k2 == i && l2 == j) {
                        swapped_term = value[n];
                        break;
                    }
		   if (i2 == l && j2 == j && k2 == i && l2 == k) {
                        swapped_term = value[n];
                        break;
                    }
		   if (i2 == j && j2 == i && k2 == k && l2 == l) {
                        swapped_term = value[n];
                        break;
                    }
		   if (i2 == k && j2 == i && k2 == j && l2 == l) {
                        swapped_term = value[n];
                        break;
                    }
		   if (i2 == k && j2 == l && k2 == j && l2 == i) {
                        swapped_term = value[n];
                        break;
                    }
		   if (i2 == j && j2 == l && k2 == k && l2 == i) {
                        swapped_term = value[n];
                        break;
                    }

                }
		//printf("swapped term= %f\n",swapped_term);

                // Apply the MP2 formula
		if (occupied[0] == occupied[1] && virtual[0] == virtual[1]) {
		    mp2_energy += integral * (2 * integral - swapped_term) / denominator; // no permutations
		} else {
		    mp2_energy += 2 * integral * (2 * integral - swapped_term) / denominator; // no permutations
		}


                // Debugging output
                //printf("MP2 Contribution: %f\n", mp2_energy);
            } else {
                // Skip terms where denominator is zero
                printf("Skipping zero denominator for (%d,%d|%d,%d), denominator = %f\n", i, j, k, l, denominator);
            }
        }
    }

    return mp2_energy;
}
