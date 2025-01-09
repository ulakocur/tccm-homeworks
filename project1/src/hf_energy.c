// hf_energy.c

#include <stdio.h>
#include "hf_energy.h"

// Function to calculate the Hartree-Fock energy
double calculate_hartree_fock_energy(double nuc_rep_energy, double* one_e_int_core, 
                                     int32_t n_up, int32_t n_orb, 
                                     int* index, double* value, int n_two_elec_int) {
    double total_energy = nuc_rep_energy;
    printf("Nucleus repulsion energy: %f\n", nuc_rep_energy);

    // One-electron energy
    double one_e_energy = 0.0;
    for (int i = 0; i < n_up; i++) {
        one_e_energy += 2.0 * one_e_int_core[i * n_orb + i]; // Diagonal elements represent <i|h|i> for one-electron integrals
    }
    total_energy += one_e_energy;
    printf("One-Electron Energy: %f\n", one_e_energy);

    // Two-electron energy
    double two_e_energy = evaluate_integrals(index, value, n_two_elec_int, n_up);
    total_energy += two_e_energy;
    printf("Two-Electron Energy: %f\n", two_e_energy);

    return total_energy;
}

// Function to evaluate two-electron integrals
double evaluate_integrals(int* index, double* value, int n_two_elec_int, int n_up) {
    double two_e_energy = 0.0;

    // Iterate over all two-electron integrals
    for (int m = 0; m < n_two_elec_int; m++) {
        int i = index[4 * m + 0];
        int j = index[4 * m + 1];
        int k = index[4 * m + 2];
        int l = index[4 * m + 3];
        double integral = value[m];

        // Only consider integrals involving occupied orbitals (spin-up electrons)
        if (i < n_up && j < n_up && k < n_up && l < n_up) {
            // Direct term: (ij|ij) - contributes positively
            if (i == k && j == l) { 
                if (i == j) { // Single permutation, not counted double
                    two_e_energy += 2 * integral;
                }
                else {
                    two_e_energy += 2 * 2 * integral; // Direct contribution (doubled to account for valid permutations)
                }
            }
            // Exchange term: (ij|ji) - contributes negatively
            if (i == j && l == k) {
                if (i == k) { // Single permutation, not counted double
                    two_e_energy -= integral;
                }
                else {
                    two_e_energy -= 2 * integral; // Exchange contribution (doubled to account for valid permutations)
                }
            }
        }
    }

    return two_e_energy;
}
