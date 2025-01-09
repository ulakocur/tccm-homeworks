# Compiler and flags
CC = gcc
CFLAGS = -O2 -Wall

# Sources and executable
SRC = main.c hf_energy.c data_gathering.c
EXTRA_SRC = mp2_hashmap.c
LOW_MEMORY_SRC = mp2_energy.c
EXEC = hf_mp2_energy.exe

# Default target
all: default

# Default version (faster)
default: $(SRC) $(EXTRA_SRC)
	$(CC) $(CFLAGS) -o $(EXEC) $(SRC) $(EXTRA_SRC) -ltrexio

# low memory version (longer computation time)
low_memory: $(SRC) $(MEM_EFFICIENT_SRC)
	$(CC) $(CFLAGS) -o $(EXEC) $(SRC) $(LOW_MEMORY_SRC) -ltrexio

# Clean up compiled files
clean:
	rm -f $(EXEC)

