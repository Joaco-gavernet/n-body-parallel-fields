#include "utils.h"

// Time measurement function
double dwalltime() {
    double sec;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

// Global variables for torus initialization
static double toroide_alfa, toroide_theta, toroide_incremento, toroide_lado, toroide_r, toroide_R;

// Memory management functions
void allocateMemory(cuerpo_t **cuerpos, float ***fuerzasX, float ***fuerzasY, float ***fuerzasZ, int N, int num_threads) {
    *cuerpos = (cuerpo_t *)malloc(N * sizeof(cuerpo_t));
    *fuerzasX = (float **)malloc(num_threads * sizeof(float *));
    *fuerzasY = (float **)malloc(num_threads * sizeof(float *));
    *fuerzasZ = (float **)malloc(num_threads * sizeof(float *));
    
    for (int i = 0; i < num_threads; i++) {
        (*fuerzasX)[i] = (float *)malloc(N * sizeof(float));
        (*fuerzasY)[i] = (float *)malloc(N * sizeof(float));
        (*fuerzasZ)[i] = (float *)malloc(N * sizeof(float));
    }
}

void freeMemory(cuerpo_t *cuerpos, float ***fuerzasX, float ***fuerzasY, float ***fuerzasZ, int num_threads) {
    free(cuerpos);
    for (int i = 0; i < num_threads; i++) {
        free((*fuerzasX)[i]);
        free((*fuerzasY)[i]);
        free((*fuerzasZ)[i]);
    }
    free(*fuerzasX);
    free(*fuerzasY);
    free(*fuerzasZ);
}

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

void inicializarCuerpos(cuerpo_t *cuerpos, int N, int num_threads, float **fuerzasX, float **fuerzasY, float **fuerzasZ) {
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

