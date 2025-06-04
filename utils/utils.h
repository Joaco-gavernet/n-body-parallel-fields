#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

// Constants for gravitational algorithm
#define PI (3.141592653589793L)
#define G 6.673e-11L
#define ESTRELLA 0
#define POLVO 1
#define H2 2

typedef struct {
    double x, y, z;
} T3Dpoint;

typedef long double tipo_ops;

// Structure representing a body
typedef struct cuerpo {
    tipo_ops masa;
    tipo_ops px;
    tipo_ops py;
    tipo_ops pz;
    tipo_ops vx;
    tipo_ops vy;
    tipo_ops vz;
    tipo_ops r;
    tipo_ops g;
    tipo_ops b;
    int tipo;
} cuerpo_t;

// Time measurement function
double dwalltime();

// Initialization functions
void inicializarEstrella(cuerpo_t *c, int i, double n);
void inicializarPolvo(cuerpo_t *c, int i, double n);
void inicializarH2(cuerpo_t *c, int i, double n);
void inicializarCuerpos(cuerpo_t *cuerpos, int N, int P); 
void inicializarFuerzas(int N, int P, tipo_ops **fuerzasX, tipo_ops **fuerzasY, tipo_ops **fuerzasZ); 


// Memory management functions
void allocateMemory(cuerpo_t **cuerpos, tipo_ops ***fuerzasX, tipo_ops ***fuerzasY, tipo_ops ***fuerzasZ, int N, int num_threads);
void freeMemory(cuerpo_t *cuerpos, tipo_ops ***fuerzasX, tipo_ops ***fuerzasY, tipo_ops ***fuerzasZ, int num_threads); 

#endif // UTILS_H
