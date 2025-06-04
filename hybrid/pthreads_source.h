#ifndef __PTHREADS_FUNCTIONS__
#define __PTHREADS_FUNCTIONS__
#include "../utils/utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>

void pthread_function(int _rank, int _N, cuerpo_t *_cuerpos, int _T, int _delta_tiempo, int _pasos); 

#endif