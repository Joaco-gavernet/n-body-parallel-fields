#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>
#include <mpi.h>
#include "../utils/utils.h"
#include "pthreads_source.h"

double *matriz_fuerzaX_local, *matriz_fuerzaY_local, *matriz_fuerzaZ_local;
double *fuerzaX_total, *fuerzaY_total, *fuerzaZ_total;
double *fuerzaX_externa, *fuerzaY_externa, *fuerzaZ_externa;

// TODO: check if this variables should be extern or not
cuerpo_t *cuerpos;
int N;
int rank;
int block_size;
int P;
int pasos;
int delta_tiempo;
pthread_barrier_t barrier;

void calcularFuerzas(cuerpo_t *cuerpos, int N, int dt, int t_id) {
  int cuerpo1, cuerpo2;
  float dif_X, dif_Y, dif_Z;
  float distancia;
  double F;
  int hasta_cuerpo = (rank == 0) ? (block_size) : N - 1;

  for (cuerpo1 = rank * (N - block_size) + t_id; cuerpo1 < hasta_cuerpo; cuerpo1 += P) {
    int id_c1 = t_id * N + cuerpo1;
    for (cuerpo2 = cuerpo1 + 1; cuerpo2 < N; cuerpo2++) {
      if ((cuerpos[cuerpo1].px == cuerpos[cuerpo2].px) && (cuerpos[cuerpo1].py == cuerpos[cuerpo2].py) && (cuerpos[cuerpo1].pz == cuerpos[cuerpo2].pz))
        continue;

      int id_c2 = t_id * N + cuerpo2;

      dif_X = cuerpos[cuerpo2].px - cuerpos[cuerpo1].px;
      dif_Y = cuerpos[cuerpo2].py - cuerpos[cuerpo1].py;
      dif_Z = cuerpos[cuerpo2].pz - cuerpos[cuerpo1].pz;

      distancia = sqrt(dif_X * dif_X + dif_Y * dif_Y + dif_Z * dif_Z);

      F = (G * cuerpos[cuerpo1].masa * cuerpos[cuerpo2].masa) / (distancia * distancia);

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

void moverCuerpos(cuerpo_t *cuerpos, int N, int dt, int t_id) {
  int hasta_cuerpo = (rank == 0) ? (block_size) : N;
  for (int i = rank * (N - block_size) + t_id; i < hasta_cuerpo; i += P) {
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

void sumarFuerzasTotales(int t_id) {
  for (int i = rank * (N - block_size) + t_id; i < N; i += P) {
    fuerzaX_total[i] = 0.0;
    fuerzaY_total[i] = 0.0;
    fuerzaZ_total[i] = 0.0;

    for (int j = 0; j < P; j++) {
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

void *thread(void *args) {
  int id = *(int *)args;

  for (int paso = 0; paso < pasos; paso++) {
    if (id == 0) {
      if (rank == 0) {
        MPI_Recv(&cuerpos[block_size], (N - block_size) * sizeof(cuerpo_t), MPI_BYTE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      } else {
        MPI_Send(&cuerpos[N - block_size], block_size * sizeof(cuerpo_t), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
      }
    }

    pthread_barrier_wait(&barrier);
    calcularFuerzas(cuerpos, N, delta_tiempo, id);
    pthread_barrier_wait(&barrier);

    sumarFuerzasTotales(id);
    pthread_barrier_wait(&barrier);
    if (rank == 0) {
      if (id == 0) {
        MPI_Send(&fuerzaX_total[block_size], (N - block_size), MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
        MPI_Send(&fuerzaY_total[block_size], (N - block_size), MPI_DOUBLE, 1, 2, MPI_COMM_WORLD);
        MPI_Send(&fuerzaZ_total[block_size], (N - block_size), MPI_DOUBLE, 1, 3, MPI_COMM_WORLD);
      }
    }
    else {
      if (id == 0) {
        MPI_Recv(fuerzaX_externa, block_size, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(fuerzaY_externa, block_size, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(fuerzaZ_externa, block_size, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      pthread_barrier_wait(&barrier);

      for (int i = (N - block_size) + id; i < N; i += P) {
        fuerzaX_total[i] += fuerzaX_externa[i - (N - block_size)];
        fuerzaY_total[i] += fuerzaY_externa[i - (N - block_size)];
        fuerzaZ_total[i] += fuerzaZ_externa[i - (N - block_size)];
      }
    }

    pthread_barrier_wait(&barrier);
    moverCuerpos(cuerpos, N, delta_tiempo, id);
    pthread_barrier_wait(&barrier);
  }

  return NULL; 
}


void finalizar() {
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

void inicializar() {
  // Initialize local force arrays
  matriz_fuerzaX_local = (double *)malloc(sizeof(double) * N * P);
  matriz_fuerzaY_local = (double *)malloc(sizeof(double) * N * P);
  matriz_fuerzaZ_local = (double *)malloc(sizeof(double) * N * P);

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
    for (int j = 0; j < P; j++) {
      matriz_fuerzaX_local[j * N + i] = 0.0;
      matriz_fuerzaY_local[j * N + i] = 0.0;
      matriz_fuerzaZ_local[j * N + i] = 0.0;
    }
  }
}

void pthread_worker(int _rank, int _N, cuerpo_t *_cuerpos, int _P, int _delta_tiempo, int _pasos) {
  N = _N;
  rank = _rank;
  // printf("From pthread_worker rank = %d\n", rank);

  if (rank == 0) block_size = 0.25 * N;
  else block_size = 0.75 * N;

  P = _P;
  pthread_t threads[P];
  int ids[P];
  pasos = _pasos;
  delta_tiempo = _delta_tiempo;
  cuerpos = _cuerpos;
  pthread_barrier_init(&barrier, NULL, P);
  inicializar(); 

  for (int i = 0; i < P; i++) {
    ids[i] = i;
    pthread_create(&threads[i], NULL, &thread, &ids[i]);
  }

  for (int i = 0; i < P; i++) pthread_join(threads[i], NULL);

  if (rank == 1) {
    MPI_Recv(cuerpos, (N - block_size) * sizeof(cuerpo_t), MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Send(cuerpos, block_size * sizeof(cuerpo_t), MPI_BYTE, 1, 0, MPI_COMM_WORLD);
  }
  finalizar(); 
}