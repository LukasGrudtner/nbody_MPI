#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

/*
 * pRNG based on http://www.cs.wm.edu/~va/software/park/park.html
 */
#define MODULUS    2147483647
#define MULTIPLIER 48271
#define DEFAULT    123456789

// sem_t maxF, cheio, newPos;
static long seed = DEFAULT;
double dt, dt_old;        /* Alterado de static para global */
double max_f;
double sim_t;             /* Simulation time */
int npart, contador;
int cnt;                  /* number of times in loop */
int tamanho_laco;         /* tamanho dos blocos de laço dos processos*/
int numt;                 /* numero de processos */
int rank;

double Random(void)
/* ----------------------------------------------------------------
 * Random returns a pseudo-random real number uniformly distributed
 * between 0.0 and 1.0.
 * ----------------------------------------------------------------
 */
{
  const long Q = MODULUS / MULTIPLIER;
  const long R = MODULUS % MULTIPLIER;
        long t;

  t = MULTIPLIER * (seed % Q) - R * (seed / Q);
  if (t > 0)
    seed = t;
  else
    seed = t + MODULUS;
  return ((double) seed / MODULUS);
}

/*
 * End of the pRNG algorithm
 */

typedef struct {
    double x, y, z;
    double mass;
} Particle;

typedef struct {
    double xold, yold, zold;
    double fx, fy, fz;
} ParticleV;

/*! para poder ser acessado no metodo computeForces
* tem que ficar localizado abaixo da struct
*/
Particle  * particles;   /* Particles */
ParticleV * pv;          /* Particle velocity */
//***************************************************

void ComputeForces()
{
  int loop = cnt;

  int i = rank * tamanho_laco;
  int backup = i;
  int aux = i + tamanho_laco;
  int sobra = npart%numt;
  int reexecutar_loop = 0;

  particles = (Particle *) malloc(sizeof(Particle)*npart);
  if((rank) == numt-1){
    aux += npart%numt;
    int tamanho_memoria_ultimo_processo = tamanho_laco + (npart%numt);
    pv = (ParticleV *) malloc(sizeof(ParticleV)*tamanho_memoria_ultimo_processo);
    MPI_Recv(&pv, tamanho_memoria_ultimo_processo, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    pv = (ParticleV *) malloc(sizeof(ParticleV)*tamanho_laco);
    MPI_Recv(&pv, tamanho_laco, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  MPI_Recv(&particles, tamanho_laco, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  Particle * others = particles;
  Particle * myparticles = particles;

  while(loop--){
    i = backup;

    for(i; i < aux; i++){
      int j;
      double xi, yi, mi, rx, ry, mj, r, fx, fy, rmin;
      rmin = 100.0;
      xi   = myparticles[i].x;
      yi   = myparticles[i].y;
      fx   = 0.0;
      fy   = 0.0;

      for (j = 0; j < npart; j++) {
        rx = xi - others[j].x;
        ry = yi - others[j].y;
        mj = others[j].mass;
        r  = rx * rx + ry * ry;
        /* ignore overlap and same particle */
        if (r == 0.0) continue;
        if (r < rmin) rmin = r;
        r  = r * sqrt(r);
        fx -= mj * rx / r;
        fy -= mj * ry / r;
      }

      pv[i].fx += fx;
      pv[i].fy += fy;
      fx = sqrt(fx*fx + fy*fy)/rmin;

      if (fx > max_f)
        max_f = fx;
    }

    MPI_Send(&max_f, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
    if (rank == numt-1) {
        MPI_Send(&pv, tamanho_laco+sobra, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
        MPI_Recv(&pv, tamanho_laco+sobra, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Send(&pv, tamanho_laco, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
        MPI_Recv(&pv, tamanho_laco, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    MPI_Recv(&particles, tamanho_laco, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Aguarda confirmação do mestre para reexecutar
    MPI_Recv(&reexecutar_loop, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}

void ComputeNewPos()
{
  int loop = cnt;
  int sobra = npart%numt;
  double max_f = 0, max_f_auxiliar = 0;
  int reexecutar_loop = 1;

  while(loop--) {

    // Recebe o pv calculado de todos os processos
    for (int rank = 1; rank < npart-1; rank++) {
        MPI_Recv(&pv[rank*tamanho_laco], sizeof(ParticleV)*tamanho_laco, MPI_CHAR, rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&max_f_auxiliar, 1, MPI_DOUBLE, rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (max_f_auxiliar > max_f)
            max_f = max_f_auxiliar;
    }
    MPI_Recv(&pv[rank*tamanho_laco], sizeof(ParticleV)*(tamanho_laco+sobra), MPI_CHAR, numt-1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&max_f_auxiliar, 1, MPI_DOUBLE, numt-1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (max_f_auxiliar > max_f)
        max_f = max_f_auxiliar;

    int i;
    double a0, a1, a2;
    double dt_new;
    a0   = 2.0 / (dt * (dt + dt_old));
    a2   = 2.0 / (dt_old * (dt + dt_old));
    a1   = -(a0 + a2);

    for (i = 0; i < npart; i++) {
        double xi, yi;
        xi             = particles[i].x;
        yi             = particles[i].y;
        particles[i].x = (pv[i].fx - a1 * xi - a2 * pv[i].xold) / a0;
        particles[i].y = (pv[i].fy - a1 * yi - a2 * pv[i].yold) / a0;
        pv[i].xold     = xi;
        pv[i].yold     = yi;
        pv[i].fx       = 0;
        pv[i].fy       = 0;
    }

    dt_new = 1.0/sqrt(max_f);
    /* Set a minimum: */
    if (dt_new < 1.0e-6)
        dt_new = 1.0e-6;
    /* Modify time step */
    if (dt_new < dt) {
        dt_old = dt;
        dt     = dt_new;
    } else if (dt_new > 4.0 * dt) {
        dt_old = dt;
        dt    *= 2.0;
    }

    sim_t += dt_old;

    for (int rank = 1; rank < numt-1; rank++) {
        MPI_Send(&particles, sizeof(Particle)*npart, MPI_CHAR, rank, MPI_ANY_TAG, MPI_COMM_WORLD);
        MPI_Send(&pv[rank*tamanho_laco], sizeof(ParticleV)*tamanho_laco, MPI_CHAR, rank, MPI_ANY_TAG, MPI_COMM_WORLD);
    }
    MPI_Send(&particles, sizeof(Particle)*npart, MPI_CHAR, numt-1, MPI_ANY_TAG, MPI_COMM_WORLD);
    MPI_Send(&pv[rank*tamanho_laco+sobra], sizeof(ParticleV)*(tamanho_laco+sobra), MPI_CHAR, numt-1, MPI_ANY_TAG, MPI_COMM_WORLD);
  }
}

void InitParticles( Particle[], ParticleV [], int );

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int i, j;
    int tmp;
    if(argc != 3){
    	printf("Wrong number of parameters.\nUsage: nbody num_bodies timesteps\n");
    	exit(1);
  	}
  	npart = atoi(argv[1]);
  	cnt = atoi(argv[2]);
  	MPI_Comm_size(MPI_COMM_WORLD, &numt);

  	dt = 0.001;
  	dt_old = 0.001;

  	if (rank == 0) {
    	/* Allocate memory for particles */
    	particles = (Particle *) malloc(sizeof(Particle)*npart);
    	pv = (ParticleV *) malloc(sizeof(ParticleV)*npart);
	}

    int aux;
    if (npart < numt){
      aux = npart;
      tamanho_laco = 1;
    } else {
      aux = numt;
      tamanho_laco = (int)(npart / (numt-1));
    }
    int sobra = npart%numt;

    /* Generate the initial values */
    if (rank == 0) {

    	InitParticles(particles, pv, npart);
        for(int rank = 1; rank < numt-1; rank++) {
          printf("1");
            MPI_Send(&pv[rank*tamanho_laco], sizeof(ParticleV)*tamanho_laco, MPI_CHAR, rank, MPI_ANY_TAG, MPI_COMM_WORLD);
            printf("2");
            MPI_Send(&particles[rank*tamanho_laco], sizeof(Particle)*tamanho_laco, MPI_CHAR, rank, MPI_ANY_TAG, MPI_COMM_WORLD);
            printf("3");
        }
        MPI_Send(&pv[rank*tamanho_laco], sizeof(ParticleV)*(tamanho_laco+sobra), MPI_CHAR, numt-1, MPI_ANY_TAG, MPI_COMM_WORLD);
        MPI_Send(&particles[rank*tamanho_laco], sizeof(Particle)*tamanho_laco, MPI_CHAR, numt-1, MPI_ANY_TAG, MPI_COMM_WORLD);
    }

    sim_t = 0.0;

    if (rank != 0) {
    	ComputeForces();
    } else {
    	ComputeNewPos();
    }

    MPI_Finalize();

    return 0;
}

void InitParticles( Particle particles[], ParticleV pv[], int npart )
{
    int i;
    for (i = 0; i < npart; i++) {
  		particles[i].x    = Random();
  		particles[i].y    = Random();
  		particles[i].z    = Random();
  		particles[i].mass = 1.0;
  		pv[i].xold    = particles[i].x;
  		pv[i].yold    = particles[i].y;
  		pv[i].zold    = particles[i].z;
  		pv[i].fx    = 0;
  		pv[i].fy    = 0;
  		pv[i].fz    = 0;
    }
}
