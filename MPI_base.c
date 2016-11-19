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
  int aux = tamanho_laco;
  int sobra = npart%(numt-1);

  particles = (Particle *) malloc(sizeof(Particle)*npart);

  if((rank) == numt-1){
    aux += sobra;
    pv = (ParticleV *) malloc(sizeof(ParticleV)*(tamanho_laco+sobra));
    MPI_Recv(&pv[0], sizeof(ParticleV)*(tamanho_laco+sobra), MPI_CHAR, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    pv = (ParticleV *) malloc(sizeof(ParticleV)*tamanho_laco);
    MPI_Recv(&pv[0], sizeof(ParticleV)*tamanho_laco, MPI_CHAR, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  MPI_Recv(&particles[0], sizeof(Particle)*npart, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  Particle * others = particles;
  Particle * myparticles = particles;

  while(loop--){
     // printf("ComputeForces[%d] - Iniciando execução %d\n", rank, loop);
    for(int i = 0; i < aux; i++){
      int j;
      double xi, yi, mi, rx, ry, mj, r, fx, fy, rmin;
      rmin = 100.0;
      xi   = myparticles[i + ((rank-1)*tamanho_laco)].x;
      yi   = myparticles[i + ((rank-1)*tamanho_laco)].y;
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

    //printf("ComputeForces[%d] - Enviando e aguardando resposta\n\n", rank);
    if (rank == numt-1) {
        MPI_Send(&max_f, 1, MPI_DOUBLE, 0, (numt-1)+10, MPI_COMM_WORLD);
        MPI_Send(&pv[0], sizeof(ParticleV)*(tamanho_laco+sobra), MPI_CHAR, 0, (numt-1), MPI_COMM_WORLD);
        MPI_Recv(&pv[0], sizeof(ParticleV)*(tamanho_laco+sobra), MPI_CHAR, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Send(&max_f, 1, MPI_DOUBLE, 0, rank+10, MPI_COMM_WORLD);
        MPI_Send(&pv[0], sizeof(ParticleV)*tamanho_laco, MPI_CHAR, 0, rank, MPI_COMM_WORLD);
        MPI_Recv(&pv[0], sizeof(ParticleV)*tamanho_laco, MPI_CHAR, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    MPI_Recv(&particles[0], sizeof(Particle)*npart, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  if (rank != numt-1)
  printf("ComputeForces[%d] - FIM DA EXECUÇÃO. Início: %d, Fim: %d\n\n", rank, (rank-1)*tamanho_laco, (rank-1)*tamanho_laco+tamanho_laco-1);
  else
  printf("ComputeForces[%d] - FIM DA EXECUÇÃO. Início: %d, Fim: %d\n\n", rank, (rank-1)*tamanho_laco, (rank-1)*tamanho_laco+tamanho_laco+sobra-1);
}

void ComputeNewPos()
{
  int loop = cnt;
  int sobra = npart%(numt-1);
  double max_f = 0, max_f_auxiliar = 0;

  while(loop--) {
    // Recebe o pv calculado de todos os processos
    for (int rank = 1; rank < numt-1; rank++) {
        MPI_Recv(&pv[(rank-1)*tamanho_laco], sizeof(ParticleV)*tamanho_laco, MPI_CHAR, rank, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&max_f_auxiliar, 1, MPI_DOUBLE, rank, rank+10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (max_f_auxiliar > max_f)
            max_f = max_f_auxiliar;
    }
    //printf("Posição inicial: [(numt-2)*tamanho_laco] = [(%d-2)*%d] = %d\n\n", numt, tamanho_laco, (numt-2)*tamanho_laco);
    MPI_Recv(&pv[(numt-2)*tamanho_laco], sizeof(ParticleV)*(tamanho_laco+sobra), MPI_CHAR, numt-1, (numt-1), MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&max_f_auxiliar, 1, MPI_DOUBLE, numt-1, (numt-1)+10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //printf("ComputeNewPos[%d] - Fechando a execução %d\n", rank, loop);

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

    //printf("ComputeNewPos[%d] - Enviando e aguardando resposta\n\n", rank);
    for (int rank = 1; rank < numt-1; rank++) {
        MPI_Send(&pv[(rank-1)*tamanho_laco], sizeof(ParticleV)*tamanho_laco, MPI_CHAR, rank, 1, MPI_COMM_WORLD);
        MPI_Send(&particles[0], sizeof(Particle)*npart, MPI_CHAR, rank, 0, MPI_COMM_WORLD);
    }
    MPI_Send(&pv[(numt-2)*tamanho_laco], sizeof(ParticleV)*(tamanho_laco+sobra), MPI_CHAR, numt-1, 1, MPI_COMM_WORLD);
    MPI_Send(&particles[0], sizeof(Particle)*npart, MPI_CHAR, numt-1, 0, MPI_COMM_WORLD);
  }
  printf("ComputerNewPos[%d] - FIM DA EXECUÇÃO.\n\n", rank);
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
    int sobra = npart%(numt-1);

    /* Generate the initial values */
    if (rank == 0) {
        printf("Init Particles\n\n");

    	InitParticles(particles, pv, npart);

        for(int rank = 1; rank < numt-1; rank++) {
            MPI_Send(&pv[(rank-1)*tamanho_laco], sizeof(ParticleV)*tamanho_laco, MPI_CHAR, rank, 1, MPI_COMM_WORLD);
            MPI_Send(&particles[0], sizeof(Particle)*npart, MPI_CHAR, rank, 0, MPI_COMM_WORLD);
        }
        MPI_Send(&pv[(numt-2)*tamanho_laco], sizeof(ParticleV)*(tamanho_laco+sobra), MPI_CHAR, numt-1, 1, MPI_COMM_WORLD);
        MPI_Send(&particles[0], sizeof(Particle)*npart, MPI_CHAR, numt-1, 0, MPI_COMM_WORLD);
    }

    sim_t = 0.0;

     if (rank != 0)  {
         printf("ComputeForces[%d]\n\n", rank);
      	ComputeForces();
    } else {
        printf("ComputeNewPos[%d]\n\n", rank);
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
