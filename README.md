# nbody_MPI

## Compilação e Execução:

- <strong>nbody_sequencial</strong>
<br/>gcc -o nbody nbody_sequencial.c -lm
<br/>./nbody (número de partículas) (número de iterações) 
<br/>// Exemplo: ./nbody 10 10

- <strong>nbody_MPI</strong>
<br/>mpicc nbody_MPI.c -o nbody_MPI -lm
<br/>mpirun -np (número de processos) ./nbody_MPI (número de partículas) (número de iterações)
<br/>// Exemplo: mpirun -np 4 ./nbody_MPI 10 10
