# Installation instructions
## Requirements
  To install this program you need to have **C compiler (gcc)** and **Make** installed on your system

## 1. Clone the repository
  `git clone https://github.com/ulakocur/tccm-homeworks`

## 2. Navigate into the project source directory
  `cd tccm-homeworks/project3/src`

## 3. Compile the code
  Compile the program using `make` 
  
  Go back to the parent directory using `cd ..`
  
## 4. Run the program
  You can run the program in the directory in which the program is located by using `./md_simulation [path_to_input_file]`
  
  An example input file [inp.txt](tests/inp.txt) is provided and it can be run using `./md_simulation tests/inp.txt`
  
## Notes: Input file structure
  The program works with the following input file structure:
    
        [Number of atoms]
        [x-coordinate] [y-coordinate] [z-coordinate] [atomic mass]
        ...
        [x-coordinate] [y-coordinate] [z-coordinate] [atomic mass]

  **Important:** This program only works on argon atoms and atomic mass other than 39.948 will result in an error.
