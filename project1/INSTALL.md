
# Installation Instructions

To run this project locally, follow these steps on a bash-based system:

## 1. Clone the repository:
git clone https://github.com/ulakocur/tccm-homeworks

## 2. Install the necessary dependencies

### HDF5:
- Ubuntu:
    sudo apt install libhdf5-dev
- macOS:
    brew install hdf5

### TREXIO:
Download the TREXIO package from the following link:
- TREXIO v2.5.0 Release (https://github.com/TREX-CoE/trexio/releases/download/v2.5.0/trexio-2.5.0.tar.gz)

Then, install it using the following commands:

cd ~/Downloads
tar -zxvf trexio-2.5.0.tar.gz
cd trexio-2.5.0
./configure
make
sudo make install

## 2. Navigate into the project source directory:

cd tccm-homeworks/project1/src

## 4. Compile the program

After installation, compile the program using the following command:

gcc -o hf_mp2_energy.exe main.c hf_energy.c data_gathering.c mp2_hashmap.c -ltrexio

## 5. Run the program

Run the compiled program with the following command, replacing '/path/to/input_file.h5' with the path to your actual input file:

./hf_mp2_energy.exe '/path/to/input_file.h5'

## 6. Optional: Compile an alternative version with less memory usage but longer computation time

If you prefer a version that uses less memory but takes longer to compute, use the following command to compile:

gcc -o hf_mp2_energy.exe main.c hf_energy.c data_gathering.c mp2_energy.c -ltrexio

## Notes

- Ensure that all dependencies (HDF5 and TREXIO) are installed correctly before proceeding with the compilation.
- If you encounter any issues during installation, please refer to the relevant documentation for HDF5 and TREXIO or open an issue in this repository.
- The directory tccm-homework/project1/tests contains various simple molecules for which the program can be tested. Expected energy values are available in the README.org file within said directory.
