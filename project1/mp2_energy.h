// mp2_energy.h

#ifndef MP2_ENERGY_H
#define MP2_ENERGY_H

// Function to calculate the MP2 energy
double calculate_mp2_energy(int* index, double* value, int n_two_elec_int, 
                             int n_up, int n_mo, double* orbital_energies);

#endif // MP2_ENERGY_H

