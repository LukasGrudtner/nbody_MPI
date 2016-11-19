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

static long seed = DEFAULT;
double dt, dt_old;        /* Alterado de static para global */
double max_f;
double sim_t;             /* Simulation time */
int npart, contador;
int cnt;                  /* number of times in loop */
int tamanho_laco;         /* tamanho dos blocos de laço dos processos*/
int nump;                 /* numero de processos */
int rank;                 /* rank do processo */

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

void ComputeForces()
{
  int loop = cnt;
  int bloco = tamanho_laco;
  int sobra = npart%(nump-1);

  if((rank) == nump-1)
    bloco += sobra;

  pv = (ParticleV *) malloc(sizeof(ParticleV)*bloco);
  particles = (Particle *) malloc(sizeof(Particle)*npart);

  Particle * others = particles;
  Particle * myparticles = particles;

  // Aguarda recebimento do primeiro bloco
  MPI_Recv(&pv[0], sizeof(ParticleV)*bloco, MPI_CHAR, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(&particles[0], sizeof(Particle)*npart, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  while (loop--){

    for (int i = 0; i < bloco; i++){
      int j;
      double xi, yi, mi, rx, ry, mj, r, fx, fy, rmin;
      rmin = 100.0;
      xi   = myparticles[i + (rank-1)*tamanho_laco].x;
      yi   = myparticles[i + (rank-1)*tamanho_laco].y;
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

    // Envia o bloco calculado para o processo mestre
    MPI_Send(&max_f, 1, MPI_DOUBLE, 0, rank+10, MPI_COMM_WORLD);
    MPI_Send(&pv[0], sizeof(ParticleV)*bloco, MPI_CHAR, 0, rank, MPI_COMM_WORLD);

    // Aguarda o recebimento do novo bloco
    MPI_Recv(&pv[0], sizeof(ParticleV)*bloco, MPI_CHAR, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&particles[0], sizeof(Particle)*npart, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  printf("ComputeForces[%d] - FIM DA EXECUÇÃO. Início: %d, Fim: %d\n", rank, (rank-1)*tamanho_laco, (rank-1)*tamanho_laco+bloco-1);
}

void ComputeNewPos()
{
  int loop = cnt;
  int sobra = npart%(nump-1);
  double max_f = 0, max_f_auxiliar = 0;

  while(loop--) {
    // Recebe o pv calculado de todos os processos
    for (int rank = 1; rank < nump-1; rank++) {
        MPI_Recv(&pv[(rank-1)*tamanho_laco], sizeof(ParticleV)*tamanho_laco, MPI_CHAR, rank, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&max_f_auxiliar, 1, MPI_DOUBLE, rank, rank+10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (max_f_auxiliar > max_f)
            max_f = max_f_auxiliar;
    }
    MPI_Recv(&pv[(nump-2)*tamanho_laco], sizeof(ParticleV)*(tamanho_laco+sobra), MPI_CHAR, nump-1, nump-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&max_f_auxiliar, 1, MPI_DOUBLE, nump-1, (nump-1)+10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

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

    // Envia os novos cálculos para todos os processos servos
    for (int rank = 1; rank < nump-1; rank++) {
        MPI_Send(&pv[(rank-1)*tamanho_laco], sizeof(ParticleV)*tamanho_laco, MPI_CHAR, rank, 1, MPI_COMM_WORLD);
        MPI_Send(&particles[0], sizeof(Particle)*npart, MPI_CHAR, rank, 0, MPI_COMM_WORLD);
    }
    MPI_Send(&pv[(nump-2)*tamanho_laco], sizeof(ParticleV)*(tamanho_laco+sobra), MPI_CHAR, nump-1, 1, MPI_COMM_WORLD);
    MPI_Send(&particles[0], sizeof(Particle)*npart, MPI_CHAR, nump-1, 0, MPI_COMM_WORLD);
  }
  printf("ComputeNewPos[%d] - FIM DA EXECUÇÃO.\n", rank);
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
  	MPI_Comm_size(MPI_COMM_WORLD, &nump);

  	dt = 0.001;
  	dt_old = 0.001;

  	if (rank == 0) {
    	/* Allocate memory for particles */
    	particles = (Particle *) malloc(sizeof(Particle)*npart);
    	pv = (ParticleV *) malloc(sizeof(ParticleV)*npart);
	}

    int bloco;
    if (npart < nump){
      bloco = npart;
      tamanho_laco = 1;
    } else {
      bloco = nump;
      tamanho_laco = (int)(npart / (nump-1));
    }
    int sobra = npart%(nump-1);

    // Inicializa as partículas e envia os primeiros blocos de cálculo para os processos servos
    if (rank == 0) {
    	InitParticles(particles, pv, npart);

        for(int rank = 1; rank < nump-1; rank++) {
            MPI_Send(&pv[(rank-1)*tamanho_laco], sizeof(ParticleV)*tamanho_laco, MPI_CHAR, rank, 1, MPI_COMM_WORLD);
            MPI_Send(&particles[0], sizeof(Particle)*npart, MPI_CHAR, rank, 0, MPI_COMM_WORLD);
        }
        MPI_Send(&pv[(nump-2)*tamanho_laco], sizeof(ParticleV)*(tamanho_laco+sobra), MPI_CHAR, nump-1, 1, MPI_COMM_WORLD);
        MPI_Send(&particles[0], sizeof(Particle)*npart, MPI_CHAR, nump-1, 0, MPI_COMM_WORLD);
    }

    sim_t = 0.0;

     if (rank != 0)  {
      	ComputeForces();
    } else {
      	ComputeNewPos();
    }

    MPI_Finalize();

    // Impressão dos valores finais das partículas
    // if (rank == 0)
    // for (i=0; i<npart; i++)
    //   fprintf(stdout,"%.5lf %.5lf\n", particles[i].x, particles[i].y);

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
