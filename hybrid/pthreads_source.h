#ifndef __PTHREADS_FUNCTIONS__
#define __PTHREADS_FUNCTIONS__
#include "../utils/utils.h"

// External declarations for shared variables
extern cuerpo_t *cuerpos;
extern int N;
extern int pasos;
extern float delta_tiempo;
extern int P;
extern pthread_barrier_t barrier;
extern float **fuerzasX, **fuerzasY, **fuerzasZ;

void pthread_worker(int rank_value, int n_value, cuerpo_t *cuerpos_values, int num_threads, int delta_t, int pasos_value);

#endif