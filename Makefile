# Unified Makefile for N-body simulation
# Compiles sequential, pthread, and hybrid (MPI+pthread) versions

# Compilers and flags
CC = gcc
MPICC = mpicc
CFLAGS = -Wall -O3 -pthread
LDFLAGS = -lm

# Get OpenMPI paths
MPI_INC := $(shell mpicc --showme:incdirs)
MPI_LIB := $(shell mpicc --showme:libdirs)

# Add MPI paths to compiler flags
MPI_CFLAGS = $(addprefix -I,$(MPI_INC))
MPI_LDFLAGS = $(addprefix -L,$(MPI_LIB))

# Common object files
UTILS_OBJ = utils/utils.o

# Targets
.PHONY: all sequential pthread hybrid clean

all: sequential pthread hybrid
	@echo "Build completed (all targets: sequential, pthread, hybrid)."

# Sequential version
sequential: $(UTILS_OBJ)
	@echo "Building sequential version..."
	$(CC) $(CFLAGS) sequential/n_body.c -o sequential/n_body $(LDFLAGS)

# Pthread version
pthread: $(UTILS_OBJ)
	@echo "Building pthread version..."
	$(CC) $(CFLAGS) pthread/n_body_pthreads.c utils/utils.c -o pthread/n_body_pthreads $(LDFLAGS)

# Hybrid version (MPI + Pthreads)
hybrid: hybrid/hybrid

hybrid/hybrid: hybrid/mpi_source.c hybrid/pthreads_source.c $(UTILS_OBJ)
	@echo "Building hybrid version (MPI + Pthreads)..."
	$(MPICC) $(CFLAGS) $(MPI_CFLAGS) -Iutils -o hybrid/hybrid \
		hybrid/mpi_source.c hybrid/pthreads_source.c utils/utils.c \
		$(LDFLAGS) $(MPI_LDFLAGS)

# Utils object file
$(UTILS_OBJ): utils/utils.c
	$(CC) $(CFLAGS) -c utils/utils.c -o $(UTILS_OBJ)

# Clean all built files
clean:
	@echo "Cleaning build files..."
	rm -f sequential/n_body pthread/n_body_pthreads hybrid/hybrid $(UTILS_OBJ) 