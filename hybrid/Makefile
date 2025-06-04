CC = mpicc
CFLAGS = -Wall -O3 -pthread -I../utils
LDFLAGS = -lm

# Get OpenMPI paths
MPI_INC := $(shell mpicc --showme:incdirs)
MPI_LIB := $(shell mpicc --showme:libdirs)

# Add MPI paths to compiler flags
CFLAGS += $(addprefix -I,$(MPI_INC))
LDFLAGS += $(addprefix -L,$(MPI_LIB))

TARGET = hybrid
SRCS = mpi_source.c pthreads_source.c
OBJS = $(SRCS:.c=.o)
UTILS_OBJ = ../utils/utils.o

all: $(TARGET)

$(TARGET): $(OBJS) $(UTILS_OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

$(UTILS_OBJ): ../utils/utils.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean 