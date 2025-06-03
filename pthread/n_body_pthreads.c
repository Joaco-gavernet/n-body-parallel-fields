// n_body_pthreads.c
// Paralelización del problema de los N-Cuerpos gravitacionales usando Pthreads
// Compilar: gcc -pthread n_body_pthreads.c ../utils/utils.c -lm -o n_body_pthreads
// Ejecutar: ./n_body_pthreads <nro_de_cuerpos> <DT> <pasos> <P>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>
#include "../utils/utils.h"


// Variables globales compartidas
cuerpo_t *cuerpos;
int N;                       // número de cuerpos
float delta_tiempo;          // intervalo de tiempo (DT)
int pasos;                   // número de pasos de simulación
int P;             // número de hilos a usar
int bodies_per_thread;        // número de cuerpos por hilo
pthread_barrier_t barrier;   // barrera para sincronizar fases

float **fuerzasX, **fuerzasY, **fuerzasZ; // Matriz de fuerzas
double toroide_alfa, toroide_theta, toroide_incremento, toroide_lado, toroide_r, toroide_R;


void calculateForces(int idW) {
    int cuerpo1, cuerpo2;
	float dif_X, dif_Y, dif_Z;
	double distancia, F;

    for (int i = idW; i < N -1; i += P) {
        // Calcular fuerzas de accion-reaccion y almacenar en fuerzas[id][j]
        for (int j = i+1; j < N; j++) {
            if ((cuerpos[i].px == cuerpos[j].px) && (cuerpos[i].py == cuerpos[j].py) && (cuerpos[i].pz == cuerpos[j].pz))
                continue;
            dif_X = cuerpos[j].px - cuerpos[i].px;
            dif_Y = cuerpos[j].py - cuerpos[i].py;
            dif_Z = cuerpos[j].pz - cuerpos[i].pz;

            distancia = sqrt(dif_X*dif_X + dif_Y*dif_Y + dif_Z*dif_Z);

            F = G * cuerpos[i].masa * cuerpos[j].masa / (distancia * distancia);

            // Cuerpo i
            fuerzasX[idW][i] += dif_X * F;
            fuerzasY[idW][i] += dif_Y * F;
            fuerzasZ[idW][i] += dif_Z * F;

            // Cuerpo j
            fuerzasX[idW][j] -= dif_X * F;
            fuerzasY[idW][j] -= dif_Y * F;
            fuerzasZ[idW][j] -= dif_Z * F;
        }
    }
}

void moveBodies(int idW) {
    double force_x, force_y, force_z;
    double dv_x, dv_y, dv_z;
    double dp_x, dp_y, dp_z;

    force_x = force_y = force_z = 0.0;
    dv_x = dv_y = dv_z = 0.0;
    dp_x = dp_y = dp_z = 0.0;

    for (int i = idW; i < N -1; i += P) {
        for (int k = 0; k < P; k++) {
            force_x += fuerzasX[k][i];
            force_y += fuerzasY[k][i];
            force_z += fuerzasZ[k][i];
            fuerzasX[k][i] = fuerzasY[k][i] = fuerzasZ[k][i] = 0.0;
        }
        dv_x = (force_x / cuerpos[i].masa) * delta_tiempo;
        dv_y = (force_y / cuerpos[i].masa) * delta_tiempo;
        dv_z = (force_z / cuerpos[i].masa) * delta_tiempo;

        dp_x = (cuerpos[i].vx + dv_x/2) * delta_tiempo;
        dp_y = (cuerpos[i].vy + dv_y/2) * delta_tiempo;
        dp_z = (cuerpos[i].vz + dv_z/2) * delta_tiempo;

        cuerpos[i].vx += dv_x;
        cuerpos[i].vy += dv_y;
        cuerpos[i].vz += dv_z;

        cuerpos[i].px += dp_x;
        cuerpos[i].py += dp_y;
        cuerpos[i].pz += dp_z;

        force_x = force_y = force_z = 0.0;
    }
}

void* simulate(void *arg) {
    int id = *(int *)arg;

    for (int paso = 0; paso < pasos; paso++) {
        calculateForces(id);
        pthread_barrier_wait(&barrier); // Sincronizar: esperar a que todos los hilos terminen cálculo de fuerzas
        moveBodies(id);
        pthread_barrier_wait(&barrier); // Sincronizar: esperar a que todos los hilos terminen de mover cuerpos
    }

    pthread_exit(NULL);
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

    allocateMemory(&cuerpos, &fuerzasX, &fuerzasY, &fuerzasZ, N, P); 
    inicializarCuerpos(cuerpos, N, P);
    inicializarFuerzas(N, P, fuerzasX, fuerzasY, fuerzasZ); 

    if (pthread_barrier_init(&barrier, NULL, P) != 0) {
        fprintf(stderr, "Error al inicializar la barrera\n");
        return -1;
    }

    pthread_t threads[P];
    int thread_ids[P];
    for (int t = 0; t < P; t++) thread_ids[t] = t; // Asignar t como id del thread
    double tIni = dwalltime(); // Inicializar temporizador

    for (int t = 0; t < P; t++) {
        if (pthread_create(&threads[t], NULL, &simulate, (void*)&thread_ids[t]) != 0) {
            fprintf(stderr, "Error al crear hilo %d\n", t);
            return -1;
        }
    }

    for (int t = 0; t < P; t++) {
        pthread_join(threads[t], NULL); // CHECK: &threads[t] or threads[t]
    }

    double tFin = dwalltime();
    double tTotal = tFin - tIni;
    printf("Tiempo en segundos (paralelo con %d hilos): %f\n", P, tTotal);

    freeMemory(cuerpos, &fuerzasX, &fuerzasY, &fuerzasZ, P);
    pthread_barrier_destroy(&barrier);

    return 0;
}

