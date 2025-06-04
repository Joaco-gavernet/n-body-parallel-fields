#include <mpi.h>
#include "pthreads_source.h"


double *matriz_fuerzaX_local, *matriz_fuerzaY_local, *matriz_fuerzaZ_local;
double *fuerzaX_total, *fuerzaY_total, *fuerzaZ_total;
double *fuerzaX_externa, *fuerzaY_externa, *fuerzaZ_externa;

cuerpo_t *cuerpos;
int N;
int rank;
int block_size;
int T;
int pasos;
int delta_tiempo;
pthread_barrier_t barrier;

void moveBodies(cuerpo_t *cuerpos, int N, int dt, int idT) {
  int hasta_cuerpo = (rank == 0) ? (block_size) : N;
  for (int i = rank * (N - block_size) + idT; i < hasta_cuerpo; i += T) {
    // Se reutiliza fuerza como aceleracion
    fuerzaX_total[i] *= 1 / cuerpos[i].masa;
    fuerzaY_total[i] *= 1 / cuerpos[i].masa;
    fuerzaZ_total[i] *= 1 / cuerpos[i].masa;

    // Calculo de velocidad
    cuerpos[i].vx += fuerzaX_total[i] * dt;
    cuerpos[i].vy += fuerzaY_total[i] * dt;
    cuerpos[i].vz += fuerzaZ_total[i] * dt;

    // Calculo de la posicion
    cuerpos[i].px += cuerpos[i].vx * dt;
    cuerpos[i].py += cuerpos[i].vy * dt;
    cuerpos[i].pz += cuerpos[i].vz * dt;
  }
}

void calculateForces(cuerpo_t *cuerpos, int N, int dt, int idT) {
  float dif_X, dif_Y, dif_Z;
  float distancia;
  double F;
  int hasta_cuerpo = (rank == 0 ? block_size : N - 1);

  // Se aplica stripes a nivel cuerpo
  // Se aplica proporcional a nivel bloque
  for (int i = rank * (N - block_size) + idT; i < hasta_cuerpo; i += T) {
    int id_c1 = idT * N + i;
    for (int j = i + 1; j < N; j++) {
      if ((cuerpos[i].px == cuerpos[j].px) && (cuerpos[i].py == cuerpos[j].py) && (cuerpos[i].pz == cuerpos[j].pz))
        continue;

      int id_c2 = idT * N + j;

      dif_X = cuerpos[j].px - cuerpos[i].px;
      dif_Y = cuerpos[j].py - cuerpos[i].py;
      dif_Z = cuerpos[j].pz - cuerpos[i].pz;

      distancia = sqrt(dif_X * dif_X + dif_Y * dif_Y + dif_Z * dif_Z);

      F = (G * cuerpos[i].masa * cuerpos[j].masa) / (distancia * distancia);

      dif_X *= F;
      dif_Y *= F;
      dif_Z *= F;

      matriz_fuerzaX_local[id_c1] += dif_X;
      matriz_fuerzaY_local[id_c1] += dif_Y;
      matriz_fuerzaZ_local[id_c1] += dif_Z;

      matriz_fuerzaX_local[id_c2] -= dif_X;
      matriz_fuerzaY_local[id_c2] -= dif_Y;
      matriz_fuerzaZ_local[id_c2] -= dif_Z;
    }
  }
}

void totalForces(int idT) {
  for (int i = rank * (N - block_size) + idT; i < N; i += T) {
    fuerzaX_total[i] = 0.0;
    fuerzaY_total[i] = 0.0;
    fuerzaZ_total[i] = 0.0;

    for (int j = 0; j < T; j++) {
      int inner_idx = i + j * N;
      fuerzaX_total[i] += matriz_fuerzaX_local[inner_idx];
      fuerzaY_total[i] += matriz_fuerzaY_local[inner_idx];
      fuerzaZ_total[i] += matriz_fuerzaZ_local[inner_idx];
      matriz_fuerzaX_local[inner_idx] = 0;
      matriz_fuerzaY_local[inner_idx] = 0;
      matriz_fuerzaZ_local[inner_idx] = 0;
    }
  }
}

void freeMem() {
  // Cleanup local arrays
  free(matriz_fuerzaX_local);
  free(matriz_fuerzaY_local);
  free(matriz_fuerzaZ_local);
  free(fuerzaX_total);
  free(fuerzaY_total);
  free(fuerzaZ_total);
  if (rank == 1) {
    free(fuerzaX_externa);
    free(fuerzaY_externa);
    free(fuerzaZ_externa);
  }
  free(cuerpos); 
  pthread_barrier_destroy(&barrier);
}

void allocMem() {
  // Initialize local force arrays
  matriz_fuerzaX_local = (double *)malloc(sizeof(double) * N * T);
  matriz_fuerzaY_local = (double *)malloc(sizeof(double) * N * T);
  matriz_fuerzaZ_local = (double *)malloc(sizeof(double) * N * T);

  fuerzaX_total = (double *)malloc(sizeof(double) * N);
  fuerzaY_total = (double *)malloc(sizeof(double) * N);
  fuerzaZ_total = (double *)malloc(sizeof(double) * N);

  if (rank == 1) {
    fuerzaX_externa = (double *)malloc(sizeof(double) * block_size);
    fuerzaY_externa = (double *)malloc(sizeof(double) * block_size);
    fuerzaZ_externa = (double *)malloc(sizeof(double) * block_size);
  }

  // Initialize force arrays to zero
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < T; j++) {
      matriz_fuerzaX_local[j * N + i] = 0.0;
      matriz_fuerzaY_local[j * N + i] = 0.0;
      matriz_fuerzaZ_local[j * N + i] = 0.0;
    }
  }
}

void *funcion(void *args) {
  int idT = *(int *)args;

  for (int paso = 0; paso < pasos; paso++) {
    // Etapa 1: enviar cuerpos hacia "menores ids"
    if (idT == 0) {
      if (rank == 0) {
        MPI_Recv(&cuerpos[block_size], (N - block_size) * sizeof(cuerpo_t), MPI_BYTE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      } else {
        MPI_Send(&cuerpos[N - block_size], block_size * sizeof(cuerpo_t), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
      }
    }

    // Etapa 2: 
    pthread_barrier_wait(&barrier);
    calculateForces(cuerpos, N, delta_tiempo, idT);
    
    // Etapa 2: 
    pthread_barrier_wait(&barrier);
    totalForces(idT);

    // Etapa 3: 
    pthread_barrier_wait(&barrier);
    if (rank == 0) {
      if (idT == 0) {
        MPI_Send(&fuerzaX_total[block_size], (N - block_size), MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
        MPI_Send(&fuerzaY_total[block_size], (N - block_size), MPI_DOUBLE, 1, 2, MPI_COMM_WORLD);
        MPI_Send(&fuerzaZ_total[block_size], (N - block_size), MPI_DOUBLE, 1, 3, MPI_COMM_WORLD);
      }
    } else {
      if (idT == 0) {
        MPI_Recv(fuerzaX_externa, block_size, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(fuerzaY_externa, block_size, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(fuerzaZ_externa, block_size, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      pthread_barrier_wait(&barrier);

      for (int i = (N - block_size) + idT; i < N; i += T) {
        fuerzaX_total[i] += fuerzaX_externa[i - (N - block_size)];
        fuerzaY_total[i] += fuerzaY_externa[i - (N - block_size)];
        fuerzaZ_total[i] += fuerzaZ_externa[i - (N - block_size)];
      }
    }

    pthread_barrier_wait(&barrier);
    moveBodies(cuerpos, N, delta_tiempo, idT);
    pthread_barrier_wait(&barrier);
  }

  pthread_exit(NULL); 
}

void pthread_function(int _rank, int _N, cuerpo_t *_cuerpos, int _T, int _delta_tiempo, int _pasos) {
  rank = _rank;
  N = _N;
  cuerpos = _cuerpos;
  T = _T;
  delta_tiempo = _delta_tiempo;
  pasos = _pasos;
  // printf("From pthread_worker rank = %d\n", rank);

  // TODO: check with different distribution like 1/8 vs 7/8
  // TODO: make dynammic mapping for P > 2  
  if (rank == 0) block_size = 0.25 * N;
  else block_size = 0.75 * N;
  pthread_t threads[T];
  int ids[T];

  pthread_barrier_init(&barrier, NULL, T);
  allocMem(); 

  for (int i = 0; i < T; i++) ids[i] = i;
  for (int i = 0; i < T; i++) pthread_create(&threads[i], NULL, &funcion, &ids[i]);
  for (int i = 0; i < T; i++) pthread_join(threads[i], NULL);

  if (rank == 1) {
    MPI_Recv(cuerpos, (N - block_size) * sizeof(cuerpo_t), MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Send(cuerpos, block_size * sizeof(cuerpo_t), MPI_BYTE, 1, 0, MPI_COMM_WORLD);
  }
  freeMem(); 
}