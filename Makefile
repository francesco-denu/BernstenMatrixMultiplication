# Compiler for serial programs
CC_SEQ = gcc

# Compiler for MPI programs
CC_MPI = mpicc

# Compiler flags for both compilers
CFLAGS = -Os

# Linker flags for math library
LDFLAGS = -lm

# Serial program source files
SERIAL_SRC = gen-int-matrix.c seq_mat_mult.c

# MPI program source files
MPI_SRC = cannon.c bernsten.c

# Serial executables (compiled with gcc)
SERIAL_EXEC = gen-int-matrix seq_mat_mult 

# MPI executables (compiled with mpicc)
MPI_EXEC = cannon bernsten

# Targets
all: $(SERIAL_EXEC) $(MPI_EXEC)

# Rule to build each serial executable
$(SERIAL_EXEC): %: %.c
	$(CC_SEQ) $(CFLAGS) -o $@ $<

# Rule to build each MPI executable
$(MPI_EXEC): %: %.c
	$(CC_MPI) $(CFLAGS) -o $@ $< $(LDFLAGS)

# Clean up
clean:
	rm -f $(SERIAL_EXEC) $(MPI_EXEC) *~
