// data_gathering.c

#include <stdio.h>
#include <stdlib.h>
#include "data_gathering.h"

int gather_data(const char* file_name, 
                double* nuc_rep_energy, 
                int32_t* n_orb, 
                int64_t* n_two_elec_int, 
                int32_t* n_up, 
                double** one_e_int_core, 
                int32_t** index, 
                double** value, 
                double** orbital_energies) {

    trexio_exit_code rc_open;
    // Open the TREXIO file in read mode (HDF5 format)
    trexio_t* file = trexio_open(file_name, 'r', TREXIO_HDF5, &rc_open);
    if (rc_open != TREXIO_SUCCESS) {
        fprintf(stderr, "Error opening file '%s': %d\n", file_name, rc_open);
        return 1;
    }

    // Read nucleus repulsion energy from the file
    trexio_exit_code nuc_rep_read = trexio_read_nucleus_repulsion(file, nuc_rep_energy);
    if (nuc_rep_read != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading nucleus repulsion energy: %d\n", nuc_rep_read);
        trexio_close(file);
        return 1;
    }

    // Read the number of molecular orbitals (MOs)
    trexio_exit_code n_orb_read = trexio_read_mo_num(file, n_orb);
    if (n_orb_read != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading number of orbitals: %d\n", n_orb_read);
        trexio_close(file);
        return 1;
    }

    // Read the number of non-zero two-electron integrals
    trexio_exit_code n_two_elec_int_read = trexio_read_mo_2e_int_eri_size(file, n_two_elec_int);
    if (n_two_elec_int_read != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading number of two-electron integrals: %d\n", n_two_elec_int_read);
        trexio_close(file);
        return 1;
    }

    // Read the number of spin-up electrons
    trexio_exit_code n_up_read = trexio_read_electron_up_num(file, n_up);
    if (n_up_read != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading number of spin-up electrons: %d\n", n_up_read);
        trexio_close(file);
        return 1;
    }

    // Allocate memory for one-electron integrals (core Hamiltonian)
    *one_e_int_core = malloc(*n_orb * *n_orb * sizeof(double));
    if (*one_e_int_core == NULL) {
        fprintf(stderr, "Memory allocation failed for one-electron integrals.\n");
        trexio_close(file);
        return 1;
    }
    
    // Read the one-electron integrals from the file
    trexio_exit_code one_e_int_core_read = trexio_read_mo_1e_int_core_hamiltonian(file, *one_e_int_core);
    if (one_e_int_core_read != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading one-electron integrals.\n");
        free(*one_e_int_core);
        trexio_close(file);
        return 1;
    }

    // Allocate memory for the indices of the two-electron integrals
    *index = malloc(4 * *n_two_elec_int * sizeof(int32_t));
    if (*index == NULL) {
        fprintf(stderr, "Memory allocation failed for indices.\n");
        free(*one_e_int_core);
        trexio_close(file);
        return 1;
    }

    // Allocate memory for the values of the two-electron integrals
    *value = malloc(*n_two_elec_int * sizeof(double));
    if (*value == NULL) {
        fprintf(stderr, "Memory allocation failed for values.\n");
        free(*index);
        free(*one_e_int_core);
        trexio_close(file);
        return 1;
    }

    // Read the two-electron integrals (indices and values)
    trexio_exit_code two_e_int_read = trexio_read_mo_2e_int_eri(file, 0, n_two_elec_int, *index, *value);
    if (two_e_int_read != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading two-electron integrals.\n");
        free(*index);
        free(*value);
        free(*one_e_int_core);
        trexio_close(file);
        return 1;
    }

    // Allocate memory for orbital energies
    *orbital_energies = malloc(*n_orb * sizeof(double));
    if (*orbital_energies == NULL) {
        fprintf(stderr, "Memory allocation failed for orbital energies.\n");
        free(*index);
        free(*value);
        free(*one_e_int_core);
        trexio_close(file);
        return 1;
    }

    // Read the orbital energies from the file
    trexio_exit_code orbital_energies_read = trexio_read_mo_energy(file, *orbital_energies);
    if (orbital_energies_read != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading orbital energies.\n");
        free(*orbital_energies);
        free(*index);
        free(*value);
        free(*one_e_int_core);
        trexio_close(file);
        return 1;
    }

    // Close the TREXIO file after all data has been read
    trexio_close(file);
    return 0;
}
