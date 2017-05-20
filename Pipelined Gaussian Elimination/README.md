## Parallel Processing - SOlution of Linear System of Equations using pipelined gaussian elimination

##### Commands to compile and run the code:

##### OpenMp
###### pgcc -mp -O4 OpenMp.c
###### ./a.out

##### OpenAcc
###### pgcc -fast -Minfo=accel -acc -ta=tesla,cuda8.0 OpenAcc.c
###### ./a.out

##### MPI
###### module load mpi/pgi_openmpi
###### mpicc MPI.c
###### mpicc -np 4 MPI


##### MPI + openMp
###### export OMP_THREADS= <no. of threads>
###### module load mpi/pgi_openmpi
###### mpicc -openmp -O4 MPI.c
###### mpicc -np 8 MPI
