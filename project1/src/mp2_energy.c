// mp2_energy.c

#include <stdio.h>
#include "mp2_energy.h"

// Function to calculate the MP2 energy
// Note: This is the inefficient method using a double for-loop which scales with O(n^2). 
// A more optimized method is available using a hashmap in the mp2_hashmap.c module, 
// which scales with O(1) and is more computationally efficient for large systems (although mp2_hashmap.c still scales O(n) due to a for loop over all indices.

double calculate_mp2_energy(int* index, double* value, int n_two_elec_int, 
                             int n_up, int n_mo, double* orbital_energies) {
    double mp2_energy = 0.0;

    // Iterate over all two-electron integrals
    for (int m = 0; m < n_two_elec_int; m++) {
        int i = index[4 * m + 0];  // Virtual (b, a)
        int j = index[4 * m + 1];
        int k = index[4 * m + 2];  // Occupied (j, i)
        int l = index[4 * m + 3];
        double integral = value[m];

        // MP2 only considers (ij|ab) integrals where i, j are occupied, and a, b are virtual
        if (((i >= n_up) + (j >= n_up) + (k >= n_up) + (l >= n_up)) == 2) {
            int occupied[2], virtual[2]; // To store the indices of occupied and virtual orbitals
            int occ_count = 0, virt_count = 0;
	    
            // Group all indices into a single array to determine occupied vs virtual orbitals
            int indices[4] = {i, j, k, l}; 
            for (int idx = 0; idx < 4; idx++) {
                if (indices[idx] < n_up) {
                    occupied[occ_count++] = indices[idx]; // Add to occupied orbitals
                } else {
                    virtual[virt_count++] = indices[idx]; // Add to virtual orbitals
                }
            }

            // Skip terms if occupied and virtual orbitals are not distinct
            if ((occupied[0] == i && occupied[1] == k) || 
                (occupied[0] == k && occupied[1] == i) || 
                (occupied[0] == j && occupied[1] == l) || 
                (occupied[0] == l && occupied[1] == j)) {
                continue;
            }

            // Calculate the denominator for the MP2 formula
            double denominator = orbital_energies[occupied[0]] + orbital_energies[occupied[1]] - 
                                 orbital_energies[virtual[0]] - orbital_energies[virtual[1]];

            // Ensure the denominator is non-zero to avoid division by zero (degenerate case)
            if (denominator != 0) {
                double swapped_term = 0.0;

                // Look for the exchange integral <ij|lk>
                for (int n = 0; n < n_two_elec_int; n++) {
                    int i2 = index[4 * n + 0];
                    int j2 = index[4 * n + 1];
                    int k2 = index[4 * n + 2];
                    int l2 = index[4 * n + 3];

                    // Check if the integral matches <ij|lk> or any of its variants
                    if ((i2 == i && j2 == j && k2 == l && l2 == k) || 
                        (i2 == i && j2 == k && k2 == l && l2 == j) ||
                        (i2 == l && j2 == k && k2 == i && l2 == j) || 
                        (i2 == l && j2 == j && k2 == i && l2 == k) ||
                        (i2 == j && j2 == i && k2 == k && l2 == l) || 
                        (i2 == k && j2 == i && k2 == j && l2 == l) || 
                        (i2 == k && j2 == l && k2 == j && l2 == i) ||
                        (i2 == j && j2 == l && k2 == k && l2 == i)) {
                        swapped_term = value[n];
                        break;
                    }
                }

                // Apply the MP2 formula
                if (occupied[0] == occupied[1] && virtual[0] == virtual[1]) {
                    mp2_energy += integral * (2 * integral - swapped_term) / denominator;
                    // Single permutation, not counted double
                } else {
                    mp2_energy += 2 * integral * (2 * integral - swapped_term) / denominator;
                    // Valid permutation, counted double (ij|kl) == (ji|lk)
                }
            } else {
                // Skip terms where denominator is zero (degenerate case)
                printf("Skipping zero denominator for (%d,%d|%d,%d), denominator = %f\n", i, j, k, l, denominator);
            }
        }
    }

    return mp2_energy;
}
