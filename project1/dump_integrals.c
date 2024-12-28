#include <stdio.h>
#include <trexio.h>
#include <stdlib.h>
#include <string.h>

// Function to evaluate all terms and permutations
double evaluate_integrals(int* index, double* value, int n_two_elec_int, int n_up);
// Function to evaluate the mp2 corrections
double mp2_energy_correction(int* index, double* value, int n_two_elec_int, int n_up);
double calculate_mp2_energy(int* index, double* value, int n_two_elec_int, int n_up, int n_mo, double* orbital_energies);

int main() {
    const char* file_name = "/Users/stef/Desktop/tccm-homeworks/project1/data/h2o.h5";
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

    // Write one-electron integrals to a file
    FILE* one_e_file = fopen("one_electron_integrals.txt", "w");
    if (one_e_file == NULL) {
        fprintf(stderr, "Error opening file for one-electron integrals.\n");
        free(one_e_int_core);
        trexio_close(file);
        return 1;
    }
    for (int i = 0; i < n_orb; i++) {
        for (int j = 0; j < n_orb; j++) {
            fprintf(one_e_file, "%d %d %.15f\n", i, j, one_e_int_core[i * n_orb + j]);
        }
    }
    fclose(one_e_file);

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

    // Write two-electron integrals to a file
    FILE* two_e_file = fopen("two_electron_integrals.txt", "w");
    if (two_e_file == NULL) {
        fprintf(stderr, "Error opening file for two-electron integrals.\n");
        free(index);
        free(value);
        free(one_e_int_core);
        trexio_close(file);
        return 1;
    }
    for (int64_t m = 0; m < n_two_elec_int; m++) {
        fprintf(two_e_file, "%d %d %d %d %.15f\n", index[4 * m + 0], index[4 * m + 1], index[4 * m + 2], index[4 * m + 3], value[m]);
    }
    fclose(two_e_file);

    // Free memory
    free(one_e_int_core);
    free(index);
    free(value);

    // Close the file
    trexio_exit_code rc_close = trexio_close(file);
    if (rc_close != TREXIO_SUCCESS) {
        fprintf(stderr, "Error closing file: %d\n", rc_close);
        return 1;
    }

    printf("TREXIO file '%s' closed successfully.\n", file_name);
    return 0;
}

