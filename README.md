# nbody_MPI

# Compilação e Execução:

- nbody_sequencial
gcc -o nbody nbody_sequencial.c -lm
./nbody (número de partículas) (número de iterações) // Exemplo: ./nbody 10 10

- nbody_MPI
mpicc nbody_MPI.c -o nbody_MPI -lm
mpirun -np (número de processos) ./nbody_MPI (número de partículas) (número de iterações)
// Exemplo: mpirun -np 4 ./nbody_MPI 10 10
