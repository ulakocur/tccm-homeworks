# Hartree-Fock and Møller–Plesset second order pertubation correction

A code in the programming language C was developed to compute the Hartree-Fock and Møller–Plesset second order pertubation correction energies from TREXIO output files.

The project is devidide into three parts.
- Data gathering from the TREXIO output files, done via the Data_Gathering.c module
- Computation of the Hartree-Fock energy, done via the hf_energy.c module
- Computation of the Møller–Plesset second order pertubation correction, done via the mp2_energy.c or mp2_hashmap.c modules depending on available memory.


