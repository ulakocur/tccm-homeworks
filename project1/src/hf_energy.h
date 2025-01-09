// hf_energy.h

#ifndef HF_ENERGY_H
#define HF_ENERGY_H

// Function to calculate the Hartree-Fock energy
double calculate_hartree_fock_energy(double nuc_rep_energy, double* one_e_int_core, 
                                     int32_t n_up, int32_t n_orb, 
                                     int* index, double* value, int n_two_elec_int);

// Function to evaluate two-electron integrals
double evaluate_integrals(int* index, double* value, int n_two_elec_int, int n_up);

#endif // HF_ENERGY_H

