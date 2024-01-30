#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
																/* Kadai_2  2023.12. 8. */
#define N 3948
#define eps 10e-12

float A[N][N], x[N], b[N];
int PIV[N];

void Make_Matrix(int n);
void Make_right(int n);
void GaussLeft( int n);
void GaussRight( int n ); 

int main(void)
{
	int i, j, n;
	double tic, toc, t, GFLOPS;

	Make_Matrix( N );
	Make_right( N );
	
	tic = omp_get_wtime();

//Solve Ax = b
								/* Gaussian elimination for the matrix A */
	GaussLeft( N );
								/* Gaussial elimination for the right-hand side vector b */
		GaussRight( N );
		toc = omp_get_wtime();


								/* if computed-b = original-x = b then your "Gauss" is correct. */
	/* printf("Vetctor x and b:\n");
	for (i=1; i<=N; i++) {
		printf(" %5.2f  --  %5.2f\n", x[i], b[i]);
	}
	printf("\n"); */

	n = N;
	GFLOPS = ((2.0/3.0)*1e-9*n*n*n)/(toc - tic);
	printf("Gaussian Elimination  Size = %3d, Time (sec) = %.2f, GFLOPS = %.2f\n", N, toc-tic, GFLOPS);
}

									//Forward Elimination for left-hand side
void GaussLeft( int n) 
{
	float valmax, swap, alfa, temp;
	int k, i, j;
	
	for (int k = 0; k < n; k++) {
		valmax = fabsf(A[k][k]);
		PIV[k] = k;
		for ( i = k+1; i < n; i++) {
			temp=fabsf(A[i][k]);
			if ( temp > valmax ) {
				valmax = temp;
				PIV[k] = i;
			}
		}

		if ( valmax > eps ) {
			if( PIV[k] != k ) {
				for ( j = k; j < n; j++ ) {
					swap = A[PIV[k]][j];
					A[PIV[k]][j] = A[k][j];
					A[k][j] = swap;
				}
			}

			for ( i=k+1; i < n; i ++ ) {
				alfa = -A[i][k]/A[k][k];
				A[i][k] = alfa;
				for ( j=k+1; j < n; j ++) {
					A[i][j] += alfa*A[k][j];
				}
			}
		}
		else {
			exit(1);
		}
	}
}

void GaussRight( int n )
{										
	float swap, sum;
	int k, i, j;
											//Forward Elimination
	for (k = 0; k < n; k++) {			
		if( PIV[k] != k ) {
			swap = b[PIV[k]];
			b[PIV[k]] = b[k];
			b[k] = swap;
		}

		for ( i = k+1; i < n; i++) {
			b[i] += A[i][k]*b[k];
		}		
	}
											//Backward Substitution
	for ( k = n-1; k >= 0; k--) {
		sum = 0.0;
		for ( j = k+1; j <=n; j++) {
			sum += + A[k][j]*b[j];
		}
		b[k] = ( b[k]-sum )/A[k][k];
	}
}

void Make_Matrix(int n)
{
	  						/* Make a matrix */
	int i, j;
	srand(3);
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			A[i][j] = (double)(rand( )%65535)/65535.0;
		}
	}	
}

void Make_right(int n)
{
 							 /* Set solution */
	int i;
	for (i=0; i<n; i++) {
		b[i] = (double)(rand( )%65535)/65535.0;
	}
}