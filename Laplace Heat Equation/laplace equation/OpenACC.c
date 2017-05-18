
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
void Init_matrix(int** mat, int n);
void Print_matrix(int** mat, int n);
void Floyd(int** mat, int n);

int main(void) {
int n;
   struct timeval start_time, stop_time, elapsed_time;  // timers

 printf("How many vertices?\n");
  scanf("%d", &n);
int *arr[n];
    for (int i=0; i<n; i++)
         arr[i] = (int *)malloc(n * sizeof(int));
  printf("the input matrix is:\n");
Init_matrix(arr, n);   
Print_matrix(arr, n);
  int num_threads=1;
   gettimeofday(&start_time,NULL);
   Floyd(arr, n);
   gettimeofday(&stop_time,NULL);
   printf("The solution is:\n");
   Print_matrix(arr, n);

   //free(arr);
   timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine

   printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);

   return 0;
}  /* main */

/*-------------------------------------------------------------------
 * Function:  Read_matrix
 * Purpose:   Read in the adjacency matrix
 * In arg:    n
 * Out arg:   mat
 */
void Init_matrix(int** mat, int n) {
   int i, j,val;
	int INFINTY=n-1;
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++)
         if (i == j)
            mat[i][j]=0;
         else {
            if ((i==j+1)|| (j==i+1)||((i==0)&&(j==n-1))||((i==n-1)&&(j==0)))
               mat[i][j]=1;
            else
               mat[i][j]= n;
         }
   }

}  /* Read_matrix */

/*-------------------------------------------------------------------
 * Function:  Print_matrix
 * Purpose:   Print the contents of the matrix
 * In args:   mat, n
 */
void Print_matrix(int**  mat, int n) {
   int i, j;

   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++)
            printf("%d ", mat[i][j]);
      printf("\n");
   }
}  /* Print_matrix */

/*-------------------------------------------------------------------
 * Function:    Floyd
 * Purpose:     Apply Floyd's algorithm to the matrix mat
 * In arg:      n
 * In/out arg:  mat:  on input, the adjacency matrix, on output
 *              lengths of the shortest paths between each pair of
 *              vertices.
 */
void Floyd(int**  mat, int n) {
   int k, i, j;
   for (k = 0; k < n; k++) {
	#pragma acc kernels
	for (i = 0; i < n; i++)
         for (j = 0; j < n; j++) {
                 mat[i][j] =fmin(mat[i][j], mat[i][k] + mat[k][j]);
         }
   }
}  /* Floyd */
