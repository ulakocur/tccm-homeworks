#include <stdio.h>
#include <stdlib.h>
#include "mp2_energy.h"

// Define a simple hashmap to store integrals by their indices (i, j, k, l).
typedef struct {
    int i, j, k, l;
    double value;
} IntegralEntry;

// Function to store an integral in the hashmap
// This function stores both the original integral (i,j|k,l) and the swapped integral (j,i|l,k) for later lookup. For some reason the looking up the integral from the original gave 0.0000 or the same as the original every time.
void store_integral(IntegralEntry* map, int* map_size, int i, int j, int k, int l, double value) {
    map[*map_size].i = i;
    map[*map_size].j = j;
    map[*map_size].k = k;
    map[*map_size].l = l;
    map[*map_size].value = value;
        
    (*map_size)++;

    // Store the swapped integral (j,i|l,k) to enable correct lookup later
    map[*map_size].i = j;
    map[*map_size].j = i;
    map[*map_size].k = l;
    map[*map_size].l = k;
    map[*map_size].value = value;

    (*map_size)++;
}

// Function to find an exchange integral in the hashmap
// Searches for the integral with the given indices (i,j|k,l) and returns its value, or 0.0 if not found.
double find_integral(IntegralEntry* map, int map_size, int i, int j, int k, int l) {
    for (int n = 0; n < map_size; n++) {
        // Check if the integral exists with matching indices (i,j|k,l)
        if (map[n].i == i && map[n].j == j && map[n].k == k && map[n].l == l) {
            return map[n].value;
        }
    }
    return 0.0;  // Return 0.0 if the integral is not found
}

// Function to calculate the MP2 energy
// This function computes the MP2 energy by iterating over two-electron integrals 
// and looking up the corresponding exchange integrals in a hashmap.
double calculate_mp2_energy(int* index, double* value, int n_two_elec_int, 
                             int n_up, int n_mo, double* orbital_energies) {
    double mp2_energy = 0.0;

    // Allocate space for the hashmap (array of integrals)
    // The hashmap stores both original and swapped integrals, so we need space for 2 * n_two_elec_int entries
    IntegralEntry* integral_map = (IntegralEntry*)malloc(2 * n_two_elec_int * sizeof(IntegralEntry));  
    int map_size = 0;

    // Step 1: Build the hashmap with only relevant integrals (those that satisfy the condition)
    for (int m = 0; m < n_two_elec_int; m++) {
        int i = index[4 * m + 0];
        int j = index[4 * m + 1];
        int k = index[4 * m + 2];
        int l = index[4 * m + 3];
        
        // Store integrals that satisfy the condition: 2 occupied and 2 virtual orbitals
        if (((i >= n_up) + (j >= n_up) + (k >= n_up) + (l >= n_up)) == 2) {
            store_integral(integral_map, &map_size, i, j, k, l, value[m]);
        }
    }

    // Step 2: Iterate over all two-electron integrals and calculate MP2 energy
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

            int indices[4] = {i, j, k, l}; // Group all indices into a single array
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

            // Ensure the denominator is non-zero to avoid division by zero
            if (denominator != 0) {
                // Look for the swapped integral (j,i|l,k) in the hashmap
                double swapped_term = find_integral(integral_map, map_size, j, i, k, l);  // Correct order for the swapped integral

                // Apply the MP2 formula
                if (occupied[0] == occupied[1] && virtual[0] == virtual[1]) {
                    mp2_energy += integral * (2 * integral - swapped_term) / denominator;
                } else {
                    mp2_energy += 2 * integral * (2 * integral - swapped_term) / denominator;
                }
            } else {
                // Skip terms where denominator is zero
                printf("Skipping zero denominator for (%d,%d|%d,%d), denominator = %f\n", i, j, k, l, denominator);
            }
        }
    }

    // Free the allocated memory for the hashmap
    free(integral_map);

    return mp2_energy;
}
