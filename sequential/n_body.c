// Compilar:
//		gcc -o n_body n_body.c -lm
// Ejecutar:
//		./n_body <nro de cuerpos> <DT> <Pasos>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

typedef long double tipo_ops;


//
// Para tiempo de ejecucion
//

double dwalltime()
{
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

double tIni, tFin, tTotal;

//
// Constantes para Algoritmo de gravitacion
//
#define PI (3.141592653589793)
#define G 6.673e-11
#define ESTRELLA 0
#define POLVO 1
#define H2 2 //Hidrogeno molecular

// ===============
// ===== CPU =====
// ===============

//
// Estructuras y variables para Algoritmo de gravitacion
//
typedef struct cuerpo cuerpo_t;
struct cuerpo{
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
	int cuerpo;
};

tipo_ops *fuerza_totalX,*fuerza_totalY, *fuerza_totalZ;
tipo_ops toroide_alfa;
tipo_ops toroide_theta;
tipo_ops toroide_incremento;
tipo_ops toroide_lado;
tipo_ops toroide_r;
tipo_ops toroide_R;

cuerpo_t *cuerpos;
int delta_tiempo = 1.0f; //Intervalo de tiempo, longitud de un paso
int pasos;
int N;

//
// Funciones para Algoritmo de gravitacion
//

void calcularFuerzas(cuerpo_t *cuerpos, int N, int dt){
	int cuerpo1, cuerpo2;
	tipo_ops dif_X, dif_Y, dif_Z;
	tipo_ops distancia;
	tipo_ops F;

	for(cuerpo1 = 0; cuerpo1<N-1 ; cuerpo1++){
		for(cuerpo2 = cuerpo1 + 1; cuerpo2<N ; cuerpo2++){
			if ( (cuerpos[cuerpo1].px == cuerpos[cuerpo2].px) && (cuerpos[cuerpo1].py == cuerpos[cuerpo2].py) && (cuerpos[cuerpo1].pz == cuerpos[cuerpo2].pz))
                continue;

			dif_X = cuerpos[cuerpo2].px - cuerpos[cuerpo1].px;
			dif_Y = cuerpos[cuerpo2].py - cuerpos[cuerpo1].py;
			dif_Z = cuerpos[cuerpo2].pz - cuerpos[cuerpo1].pz;
                
			distancia = sqrt(dif_X*dif_X + dif_Y*dif_Y + dif_Z*dif_Z);

			F = (G*cuerpos[cuerpo1].masa*cuerpos[cuerpo2].masa)/(distancia*distancia);

			dif_X *= F;
			dif_Y *= F;
			dif_Z *= F;

			fuerza_totalX[cuerpo1] += dif_X;
			fuerza_totalY[cuerpo1] += dif_Y;
			fuerza_totalZ[cuerpo1] += dif_Z;

			fuerza_totalX[cuerpo2] -= dif_X;
			fuerza_totalY[cuerpo2] -= dif_Y;
			fuerza_totalZ[cuerpo2] -= dif_Z;
		}
	}
}

void moverCuerpos(cuerpo_t *cuerpos, int N, int dt){
	int cuerpo;
	for(cuerpo = 0; cuerpo<N ; cuerpo++){
		fuerza_totalX[cuerpo] *= 1/cuerpos[cuerpo].masa;
		fuerza_totalY[cuerpo] *= 1/cuerpos[cuerpo].masa;
		fuerza_totalZ[cuerpo] *= 1/cuerpos[cuerpo].masa;

		cuerpos[cuerpo].vx += fuerza_totalX[cuerpo]*dt;
		cuerpos[cuerpo].vy += fuerza_totalY[cuerpo]*dt;
		cuerpos[cuerpo].vz += fuerza_totalZ[cuerpo]*dt;

		cuerpos[cuerpo].px += cuerpos[cuerpo].vx *dt;
		cuerpos[cuerpo].py += cuerpos[cuerpo].vy *dt;
		cuerpos[cuerpo].pz += cuerpos[cuerpo].vz *dt;

		fuerza_totalX[cuerpo] = 0.0;
		fuerza_totalY[cuerpo] = 0.0;
		fuerza_totalZ[cuerpo] = 0.0;
	}
}

void gravitacionCPU(cuerpo_t *cuerpos, int N, int dt){
	calcularFuerzas(cuerpos,N,dt);
	moverCuerpos(cuerpos,N,dt);
}

void inicializarEstrella(cuerpo_t *cuerpo,int i,double n){

    cuerpo->masa = 0.001*8;

	if ((toroide_alfa + toroide_incremento) >=2*M_PI){
		toroide_alfa = 0;
		toroide_theta += toroide_incremento;
	}else{
		toroide_alfa+=toroide_incremento;
	}

	cuerpo->px = (toroide_R + toroide_r*cos(toroide_alfa))*cos(toroide_theta);
	cuerpo->py = (toroide_R + toroide_r*cos(toroide_alfa))*sin(toroide_theta);
	cuerpo->pz = toroide_r*sin(toroide_alfa);

	cuerpo->vx = 0.0;
	cuerpo->vy = 0.0;
	cuerpo->vz = 0.0;

	cuerpo->r = 1.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->g = 1.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->b = 1.0; //(double )rand()/(RAND_MAX+1.0);
}

void inicializarPolvo(cuerpo_t *cuerpo,int i,double n){

    cuerpo->masa = 0.001*4;
	
        if ((toroide_alfa + toroide_incremento) >=2*M_PI){
            toroide_alfa = 0;
            toroide_theta += toroide_incremento;
        }else{
            toroide_alfa+=toroide_incremento;
        }

	cuerpo->px = (toroide_R + toroide_r*cos(toroide_alfa))*cos(toroide_theta);
	cuerpo->py = (toroide_R + toroide_r*cos(toroide_alfa))*sin(toroide_theta);
	cuerpo->pz = toroide_r*sin(toroide_alfa);
	
	cuerpo->vx = 0.0;
	cuerpo->vy = 0.0;
	cuerpo->vz = 0.0;
    
	cuerpo->r = 1.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->g = 0.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->b = 0.0; //(double )rand()/(RAND_MAX+1.0);
}

void inicializarH2(cuerpo_t *cuerpo,int i,double n){

    cuerpo->masa = 0.001;

	if ((toroide_alfa + toroide_incremento) >=2*M_PI){
            toroide_alfa = 0;
            toroide_theta += toroide_incremento;
	}else{
            toroide_alfa+=toroide_incremento;
	}

	cuerpo->px = (toroide_R + toroide_r*cos(toroide_alfa))*cos(toroide_theta);
	cuerpo->py = (toroide_R + toroide_r*cos(toroide_alfa))*sin(toroide_theta);
	cuerpo->pz = toroide_r*sin(toroide_alfa);

	cuerpo->vx = 0.0;
	cuerpo->vy = 0.0;
	cuerpo->vz = 0.0;

	cuerpo->r = 1.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->g = 1.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->b = 1.0; //(double )rand()/(RAND_MAX+1.0);
}

void inicializarCuerpos(cuerpo_t *cuerpos,int N){
	int cuerpo;
	double n = N;

	toroide_alfa = 0.0;
	toroide_theta = 0.0;
	toroide_lado = sqrt(N);
	toroide_incremento = 2*M_PI / toroide_lado;
	toroide_r = 1.0;
	toroide_R = 2*toroide_r;
	
	srand(1);  // Fixed seed for reproducibility

	for(cuerpo = 0; cuerpo < N; cuerpo++){
        fuerza_totalX[cuerpo] = 0.0;
		fuerza_totalY[cuerpo] = 0.0;
		fuerza_totalZ[cuerpo] = 0.0;

		cuerpos[cuerpo].cuerpo = (rand() %3);

		if (cuerpos[cuerpo].cuerpo == ESTRELLA){
			inicializarEstrella(&cuerpos[cuerpo],cuerpo,n);
		}else if (cuerpos[cuerpo].cuerpo == POLVO){
			inicializarPolvo(&cuerpos[cuerpo],cuerpo,n);
		}else if (cuerpos[cuerpo].cuerpo == H2){
			inicializarH2(&cuerpos[cuerpo],cuerpo,n);
		}
	}

	cuerpos[0].masa = 2.0e2;
	cuerpos[0].px = 0.0;
	cuerpos[0].py = 0.0;
	cuerpos[0].pz = 0.0;
	cuerpos[0].vx = -0.000001;
	cuerpos[0].vy = -0.000001;
	cuerpos[0].vz = 0.0;

	cuerpos[1].masa = 1.0e1;
	cuerpos[1].px = -1.0;
	cuerpos[1].py = 0.0;
	cuerpos[1].pz = 0.0;
	cuerpos[1].vx = 0.0;
	cuerpos[1].vy = 0.0001;
	cuerpos[1].vz = 0.0;
}

void finalizar(void){
	free(cuerpos);
	free(fuerza_totalX);
	free(fuerza_totalY);
	free(fuerza_totalZ);
}

int main(int argc, char * argv[]) {

	if (argc < 4){
		printf("Ejecutar: %s <nro. de cuerpos> <DT> <pasos>\n",argv[0]);
		return -1;
	}
	
	N = atoi(argv[1]);
	delta_tiempo = atof(argv[2]);
	pasos = atoi(argv[3]);
	
	cuerpos = (cuerpo_t*)malloc(sizeof(cuerpo_t)*N);
	fuerza_totalX = (tipo_ops*)malloc(sizeof(tipo_ops)*N);
	fuerza_totalY = (tipo_ops*)malloc(sizeof(tipo_ops)*N);
	fuerza_totalZ = (tipo_ops*)malloc(sizeof(tipo_ops)*N);

	inicializarCuerpos(cuerpos,N);
	
	tIni = dwalltime(); 
	
	int paso;
	for(paso=0; paso<pasos; paso++){
		gravitacionCPU(cuerpos,N,delta_tiempo);
	}

	tFin =	dwalltime();
	tTotal = tFin - tIni;
	
	printf("Tiempo en segundos: %Lf\n",tTotal);
	printf("Posiciones finales de los cuerpos:\n");
	for (int i = 0; i < N; i++) {
		printf("Cuerpo %d: (%.6Lf, %.6Lf, %.6Lf)\n", i, cuerpos[i].px, cuerpos[i].py, cuerpos[i].pz);
	}

	finalizar();
    return(0);

}
