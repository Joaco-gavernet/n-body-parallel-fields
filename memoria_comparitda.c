// n_body_pthreads.c
// Paralelización del problema de los N-Cuerpos gravitacionales usando Pthreads
// Compilar: gcc -pthread -o n_body_pthreads n_body_pthreads.c -lm
// Ejecutar: ./n_body_pthreads <nro_de_cuerpos> <DT> <pasos> <num_threads>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>

//
// Para medir tiempo de ejecución
//
double dwalltime() {
    double sec;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

//
// Constantes para el algoritmo de gravitación
//
#define PI (3.141592653589793)
#define G 6.673e-11f
#define ESTRELLA 0
#define POLVO 1
#define H2 2

//
// Estructura que representa un cuerpo
//
typedef struct cuerpo {
    float masa;
    float px;
    float py;
    float pz;
    float vx;
    float vy;
    float vz;
    float r;
    float g;
    float b;
    int tipo;
} cuerpo_t;

//
// Variables globales compartidas
//
cuerpo_t *cuerpos;
int N;                       // número de cuerpos
float delta_tiempo;          // intervalo de tiempo (DT)
int pasos;                   // número de pasos de simulación
int num_threads;             // número de hilos a usar
int bodies_per_thread;        // número de cuerpos por hilo
pthread_barrier_t barrier;   // barrera para sincronizar fases

// To-Do: Eliminar fuerza_totalX, fuerza_totalY, fuerza_totalZ
float *fuerza_totalX, *fuerza_totalY, *fuerza_totalZ;
float **fuerzasX, **fuerzasY, **fuerzasZ; // Matriz de fuerzas

//
// Funciones de inicialización (idénticas al secuencial)
//
double toroide_alfa, toroide_theta, toroide_incremento, toroide_lado, toroide_r, toroide_R;

void inicializarEstrella(cuerpo_t *c, int i, double n) {
    c->masa = 0.001f * 8.0f;

    if ((toroide_alfa + toroide_incremento) >= 2 * M_PI) {
        toroide_alfa = 0;
        toroide_theta += toroide_incremento;
    } else {
        toroide_alfa += toroide_incremento;
    }

    c->px = (toroide_R + toroide_r * cos(toroide_alfa)) * cos(toroide_theta);
    c->py = (toroide_R + toroide_r * cos(toroide_alfa)) * sin(toroide_theta);
    c->pz = toroide_r * sin(toroide_alfa);

    c->vx = 0.0f;
    c->vy = 0.0f;
    c->vz = 0.0f;

    c->r = 1.0f;
    c->g = 1.0f;
    c->b = 1.0f;
}

void inicializarPolvo(cuerpo_t *c, int i, double n) {
    c->masa = 0.001f * 4.0f;

    if ((toroide_alfa + toroide_incremento) >= 2 * M_PI) {
        toroide_alfa = 0;
        toroide_theta += toroide_incremento;
    } else {
        toroide_alfa += toroide_incremento;
    }

    c->px = (toroide_R + toroide_r * cos(toroide_alfa)) * cos(toroide_theta);
    c->py = (toroide_R + toroide_r * cos(toroide_alfa)) * sin(toroide_theta);
    c->pz = toroide_r * sin(toroide_alfa);

    c->vx = 0.0f;
    c->vy = 0.0f;
    c->vz = 0.0f;

    c->r = 1.0f;
    c->g = 0.0f;
    c->b = 0.0f;
}

void inicializarH2(cuerpo_t *c, int i, double n) {
    c->masa = 0.001f;

    if ((toroide_alfa + toroide_incremento) >= 2 * M_PI) {
        toroide_alfa = 0;
        toroide_theta += toroide_incremento;
    } else {
        toroide_alfa += toroide_incremento;
    }

    c->px = (toroide_R + toroide_r * cos(toroide_alfa)) * cos(toroide_theta);
    c->py = (toroide_R + toroide_r * cos(toroide_alfa)) * sin(toroide_theta);
    c->pz = toroide_r * sin(toroide_alfa);

    c->vx = 0.0f;
    c->vy = 0.0f;
    c->vz = 0.0f;

    c->r = 1.0f;
    c->g = 1.0f;
    c->b = 1.0f;
}

void inicializarCuerpos(cuerpo_t *cuerpos, int N) {
    int i;
    double n = (double)N;

    toroide_alfa = 0.0;
    toroide_theta = 0.0;
    toroide_lado = sqrt((double)N);
    toroide_incremento = 2 * M_PI / toroide_lado;
    toroide_r = 1.0;
    toroide_R = 2.0 * toroide_r;

    srand(time(NULL));

    for (i = 0; i < N; i++) {
        for (int j = 0; j < num_threads; j++) {
            fuerzasX[j][i] = 0.0f;
            fuerzasY[j][i] = 0.0f;
            fuerzasZ[j][i] = 0.0f;
        }

        cuerpos[i].tipo = rand() % 3;
        if (cuerpos[i].tipo == ESTRELLA) {
            inicializarEstrella(&cuerpos[i], i, n);
        } else if (cuerpos[i].tipo == POLVO) {
            inicializarPolvo(&cuerpos[i], i, n);
        } else {
            inicializarH2(&cuerpos[i], i, n);
        }
    }

    // Configuración manual de dos cuerpos para simetría inicial
    cuerpos[0].masa = 2.0e2f;
    cuerpos[0].px   = 0.0f;
    cuerpos[0].py   = 0.0f;
    cuerpos[0].pz   = 0.0f;
    cuerpos[0].vx   = -0.000001f;
    cuerpos[0].vy   = -0.000001f;
    cuerpos[0].vz   = 0.0f;

    cuerpos[1].masa = 1.0e1f;
    cuerpos[1].px   = -1.0f;
    cuerpos[1].py   =  0.0f;
    cuerpos[1].pz   =  0.0f;
    cuerpos[1].vx   =  0.0f;
    cuerpos[1].vy   =  0.0001f;
    cuerpos[1].vz   =  0.0f;
}

void calculateForces(int idW) {
    int cuerpo1, cuerpo2;
	float dif_X, dif_Y, dif_Z;
	double distancia, F;

    for (int i = idW; i < bodies_per_thread -1; i += num_threads) {
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

    for (int i = idW; i < bodies_per_thread -1; i += num_threads) {
        for (int k = 0; k < num_threads; k++) {
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

//
// Función que cada hilo ejecuta para simular
//
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
        printf("Uso: %s <nro_de_cuerpos> <DT> <pasos> <num_threads>\n", argv[0]);
        return -1;
    }

    // Leer parámetros
    N = atoi(argv[1]);
    delta_tiempo = atof(argv[2]); // DT
    pasos = atoi(argv[3]);
    num_threads = atoi(argv[4]);
    bodies_per_thread = N / num_threads;

    // Reservar memoria para cuerpos y arreglos de fuerza
    cuerpos = (cuerpo_t *)malloc(sizeof(cuerpo_t) * N);
    fuerzasX = (float **)malloc(sizeof(float *) * num_threads);
    fuerzasY = (float **)malloc(sizeof(float *) * num_threads);
    fuerzasZ = (float **)malloc(sizeof(float *) * num_threads);
    for (int i = 0; i < num_threads; i++) {
        fuerzasX[i] = (float *)malloc(sizeof(float) * N);
        fuerzasY[i] = (float *)malloc(sizeof(float) * N);
        fuerzasZ[i] = (float *)malloc(sizeof(float) * N);
    }
    if (!cuerpos || !fuerzasX || !fuerzasY || !fuerzasZ) {
        fprintf(stderr, "Error al asignar memoria\n");
        return -1;
    }
    inicializarCuerpos(cuerpos, N);

    // Inicializar barrera para sincronizar num_threads hilos
    if (pthread_barrier_init(&barrier, NULL, num_threads) != 0) {
        fprintf(stderr, "Error al inicializar la barrera\n");
        return -1;
    }

    /*/////////////////////////////////////////////////
                    1. DEFINE THREADS
    /////////////////////////////////////////////////*/

    pthread_t threads[num_threads];
    int thread_ids[num_threads];

    for (int t = 0; t < num_threads; t++) thread_ids[t] = t; // Asignar t como id del thread

    double tIni = dwalltime(); // Inicializar temporizador

    /*/////////////////////////////////////////////////
                    2. CREATE THREADS
    /////////////////////////////////////////////////*/

    for (int t = 0; t < num_threads; t++) {
        if (pthread_create(&threads[t], NULL, simulate, (void*)&thread_ids[t]) != 0) {
            fprintf(stderr, "Error al crear hilo %d\n", t);
            return -1;
        }
    }

    /*/////////////////////////////////////////////////
                    3. JOIN THREADS
    /////////////////////////////////////////////////*/

    // Esperar a que todos los hilos terminen
    for (int t = 0; t < num_threads; t++) {
        pthread_join(threads[t], NULL); // CHECK: &threads[t] or threads[t]
    }

    double tFin = dwalltime();
    double tTotal = tFin - tIni;
    printf("Tiempo en segundos (paralelo con %d hilos): %f\n", num_threads, tTotal);

    // Destruir barrera y liberar memoria
    // To-Do: verificar que se libera toda la memoria usada
    pthread_barrier_destroy(&barrier);
    free(cuerpos);
    for (int i = 0; i < num_threads; i++) {
        free(fuerzasX[i]);
        free(fuerzasY[i]);
        free(fuerzasZ[i]);
    }
    free(fuerzasX);
    free(fuerzasY);
    free(fuerzasZ);


    return 0;
}

