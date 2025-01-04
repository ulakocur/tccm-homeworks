// data_gathering.h

#ifndef DATA_GATHERING_H
#define DATA_GATHERING_H

#include <trexio.h>

// Function to gather data from the TREXIO file
int gather_data(const char* file_name, 
                double* nuc_rep_energy, 
                int32_t* n_orb, 
                int64_t* n_two_elec_int, 
                int32_t* n_up, 
                double** one_e_int_core, 
                int32_t** index, 
                double** value, 
                double** orbital_energies);

#endif // DATA_GATHERING_H

