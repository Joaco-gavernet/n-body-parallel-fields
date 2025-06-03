// TODO: add some comments on compilation and execution here...
// Compilar: make
// Ejecutar: mpirun -np 2 ./hybrid/nbody <nro_de_cuerpos> <DT> <pasos> <P>

#include <pthread.h>
#include <mpi.h>
#include "../utils/utils.h"
#include "pthreads_source.h"


// Estructura para datos del worker
typedef struct {
    int idW;           // Worker ID
    int P;            // Total number of workers
    int start;        // Start timestep
    int finish;       // End timestep
    double DT;        // Time step
    int blockSize;    // Size of local block
    int tempSize;     // Max size of temporary arrays
    T3Dpoint *m;      // Masses
    T3Dpoint *p;      // Positions
    T3Dpoint *v;      // Velocities
    T3Dpoint *f;      // Forces
    T3Dpoint *tp;     // Temporary positions
    T3Dpoint *tv;     // Temporary velocities
    T3Dpoint *tf;     // Temporary forces
} WorkerData;

// Variables globales compartidas
cuerpo_t *cuerpos;
int N;                       // número de cuerpos
float delta_tiempo;          // intervalo de tiempo (DT)
int pasos;                   // número de pasos de simulación
int P;                       // número de hilos a usar
int bodies_per_thread;       // número de cuerpos por hilo
int MPI_procesos; 
pthread_barrier_t barrier;   // barrera para sincronizar fases
float **fuerzasX, **fuerzasY, **fuerzasZ; // Matriz de fuerzas


static void mpi_main(int rank, int N, int T, int delta_tiempo, int pasos, cuerpo_t *cuerpos) {
    double tIni, tFin;

    printf("Hello from rank %d\n", rank);

    pthread_barrier_init(&barrier, NULL, P);  
    allocateMemory(&cuerpos, &fuerzasX, &fuerzasY, &fuerzasZ, N, P);
    if (rank == 0) inicializarCuerpos(cuerpos, N, P, fuerzasX, fuerzasY, fuerzasZ);

    tIni = dwalltime();
    MPI_Bcast(cuerpos, N * sizeof(cuerpo_t), MPI_BYTE, 0, MPI_COMM_WORLD);
    pthread_worker(rank, N, cuerpos, P, delta_tiempo, pasos);
    tFin = dwalltime();

    if (rank == 1) printf("Tiempo total de ejecucion: %.10f segundos\n", tFin - tIni);

    freeMemory(cuerpos, &fuerzasX, &fuerzasY, &fuerzasZ, P);
    pthread_barrier_destroy(&barrier);
}


int main(int argc, char *argv[]) {
    if (argc < 5) {
        printf("Uso: %s <nro_de_cuerpos> <DT> <pasos> <P>\n", argv[0]);
        return -1;
    }

    N = atoi(argv[1]);
    delta_tiempo = atof(argv[2]); // DT
    pasos = atoi(argv[3]);
    P = atoi(argv[4]);
    bodies_per_thread = N / P;

    int rank, size, provided; 
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided < MPI_THREAD_MULTIPLE) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    mpi_main(rank, N, P, delta_tiempo, pasos, cuerpos);

    MPI_Barrier(MPI_COMM_WORLD); 
    MPI_Finalize();

    return 0;
}

