// Compilar: make
// Ejecutar: mpirun --bind-to none -np 2 hybrid <nro_de_cuerpos> <DT> <pasos> <P>

#include <mpi.h>
#include "pthreads_source.h"


static void mpi_main(int rank, int N, int P, int delta_tiempo, int pasos, cuerpo_t *cuerpos) {
    double tIni, tFin;

    // printf("From mpi_main rank = %d\n", rank);
    if (rank == 0) inicializarCuerpos(cuerpos, N, P);

    tIni = dwalltime();
    MPI_Bcast(cuerpos, N * sizeof(cuerpo_t), MPI_BYTE, 0, MPI_COMM_WORLD);
    pthread_function(rank, N, cuerpos, P, delta_tiempo, pasos);
    tFin = dwalltime();

    if (rank == 0) {
        printf("Tiempo en segundos: %.10f\n", tFin - tIni);
        printf("Posiciones finales de los cuerpos:\n");
        for (int i = 0; i < N; i++) {
            printf("Cuerpo %d: px=%.6f, py=%.6f, pz=%.6f\n", i, cuerpos[i].px, cuerpos[i].py, cuerpos[i].pz);
        }
    }
}


int main(int argc, char *argv[]) {
    cuerpo_t *cuerpos; 
    int N, rank, P, delta_tiempo, pasos, size, provided; 

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided < MPI_THREAD_MULTIPLE) MPI_Abort(MPI_COMM_WORLD, 1);

    if (argc < 5) {
        printf("Uso: %s <nro_de_cuerpos> <DT> <pasos> <P>\n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, 1);
        return -1;
    }
    N = atoi(argv[1]);
    delta_tiempo = atof(argv[2]); // DT
    pasos = atoi(argv[3]);
    P = atoi(argv[4]);

    cuerpos = (cuerpo_t *)malloc(sizeof(cuerpo_t) * N); 

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // printf("From main size = %d and rank = %d\n", size, rank); 
    mpi_main(rank, N, P, delta_tiempo, pasos, cuerpos);

    MPI_Barrier(MPI_COMM_WORLD); 
    // printf("Checkpoint rank = %d\n", rank); 
    MPI_Finalize();

    return 0;
}

