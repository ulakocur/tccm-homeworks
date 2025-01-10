# Hartree-Fock and Møller–Plesset Second-Order Perturbation Correction

This project provides a C-based implementation for calculating the Hartree-Fock and Møller–Plesset second-order perturbation correction energies using data from TREXIO output files.

### Project Structure

- **[INSTALL.md](INSTALL.md):** Instructions for compiling and running the program.
- **[tests](tests):** Example input files for testing.
- **[src](src):** Source code of the program, including the `Makefile`.

The project is organized into three main components:

1. **Data Gathering:** Extracts required data from TREXIO output files, implemented in `Data_Gathering.c`.
2. **Hartree-Fock Energy Calculation:** Computes the Hartree-Fock energy, implemented in `hf_energy.c`.
3. **Møller–Plesset Second-Order Perturbation Correction:** Computes the MP2 energy using either:
   - `mp2_energy.c` (default approach).
   - `mp2_hashmap.c` (memory-optimized approach).

### Additional Files

- **[LICENSE](LICENSE):** Licensing information for the code.
- **[AUTHORS.md](AUTHORS.md):** List of contributors to the project.
