## Parallel Processing - Laplace Heat Equation

### Laplace Solver with CPU code for openMp and openACC.

GFLOPS for openMP only - 44 gflops
GFLOPS for openACC only - 52.64 gflops

### Command to run openMp
export OMP_NUM_THREADS=28
pgcc -mp O4 laplace_parallel_omp.c
./a.out

<p align="center">
 <img src="https://github.com/rachnasidana28/ParallelProcessing/blob/master/Laplace%20Heat%20Equation/laplace%20equation/images/openMp1.png" width="450"/>
  <img src="https://github.com/rachnasidana28/ParallelProcessing/blob/master/Laplace%20Heat%20Equation/laplace%20equation/images/openMp2.png" width="450"/>
  <img src="https://github.com/rachnasidana28/ParallelProcessing/blob/master/Laplace%20Heat%20Equation/laplace%20equation/images/openMp3.png" width="450"/>
</p>

### Command to run openACC
pgcc -acc -ta=tesla,cuda8.0 laplace_parallel_acc.c
./a.out
