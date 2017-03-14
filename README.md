# ParallelProcessing

--> Laplace Solver with CPU code for openMp and openACC.

GFLOPS for openMP only - 44 gflops
GFLOPS for openACC only - 52.64 gflops

#Command to run openMp
export OMP_NUM_THREADS=28
pgcc -mp O4 laplace_parallel_omp.c
./a.out

#Command to run openACC
pgcc -acc -ta=tesla,cuda8.0 laplace_parallel_acc.c
./a.out
