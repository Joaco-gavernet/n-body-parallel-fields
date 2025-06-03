#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

// Constants for gravitational algorithm
#define PI (3.141592653589793)
#define G 6.673e-11f
#define ESTRELLA 0
#define POLVO 1
#define H2 2

// Structure representing a body
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

// Time measurement function
double dwalltime();

// Initialization functions
void inicializarEstrella(cuerpo_t *c, int i, double n);
void inicializarPolvo(cuerpo_t *c, int i, double n);
void inicializarH2(cuerpo_t *c, int i, double n);
void inicializarCuerpos(cuerpo_t *cuerpos, int N, int num_threads, float **fuerzasX, float **fuerzasY, float **fuerzasZ);

// Memory management functions
void allocateMemory(cuerpo_t **cuerpos, float ***fuerzasX, float ***fuerzasY, float ***fuerzasZ, int N, int num_threads);
void freeMemory(cuerpo_t *cuerpos, float ***fuerzasX, float ***fuerzasY, float ***fuerzasZ, int num_threads); 

#endif // UTILS_H
