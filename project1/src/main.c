// main.c

#include <stdio.h>
#include <stdlib.h>
#include "data_gathering.h"
#include "hf_energy.h"
#include "mp2_energy.h"

int main(int argc, char* argv[]) {
    if (argc != 2) {  // Check if exactly 1 argument (path to TREXIO file) is provided
        fprintf(stderr, "Usage: %s <path_to_file>\n", argv[0]); 
        return 1;
    }

    const char* file_name = argv[1];  // Retrieve the path from the command-line argument

    // Define variables and their types for further calculations
    double nuc_rep_energy;
    int32_t n_orb;
    int64_t n_two_elec_int;
    int32_t n_up;
    double* one_e_int_core = NULL;
    int32_t* index = NULL;
    double* value = NULL;
    double* orbital_energies = NULL;

    // Gather data from the TREXIO file using the data_gathering.c module
    int result = gather_data(file_name, &nuc_rep_energy, &n_orb, &n_two_elec_int, &n_up, &one_e_int_core, &index, &value, &orbital_energies);
    if (result != 0) {
        return 1; // Error handling is done within the gather_data function
    }


    printf("\nStarting energy calculation...\n");
    printf("Hartree-Fock energy calculation starting...\n");

    // Perform Hartree-Fock energy calculation using the hf_energy.c module
    double total_energy = calculate_hartree_fock_energy(nuc_rep_energy, one_e_int_core, 
                                                       n_up, n_orb, 
                                                       index, value, n_two_elec_int);

    // Print Hartree-Fock energy result
    printf("Total Hartree-Fock Energy: %f\n\n", total_energy);
    printf("Finished Hartree-Fock energy calculation successfully!\n");
    printf("\nStarting Møller–Plesset second order energy correction calculation...\n");

    // Perform MP2 energy calculation using the mp2_hashmap.c module.
    // mp2_energy.c can also be used, but it utilizes a double for-loop (O(n^2)) instead of a hashmap (O(1)).
    double mp2_energy = calculate_mp2_energy(index, value, n_two_elec_int, n_up, n_orb, orbital_energies);
    total_energy += mp2_energy;
    printf("Møller–Plesset second order energy correction: %f\n", mp2_energy);
    printf("\nMøller–Plesset second order energy correction calculation finished successfully!\n\n");
    printf("Energy Calculation finished successfully!\n\n");

    // Print final total energy (Hartree-Fock + MP2)
    printf("Total Energy (Hartree-Fock + MP2): %f\n", total_energy);

    // Free dynamically allocated memory
    free(one_e_int_core);
    free(index);
    free(value);
    free(orbital_energies);

    printf("Thank you for using our program!\n\n");

    return 0;
}
