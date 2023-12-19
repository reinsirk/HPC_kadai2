#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <time.h>
#include <omp.h>
/* Kadai_2  2023.12. 8. */
#define N 6000

double A[N + 1][N + 1], x[N + 1], b[N + 1];
//double At[N + 1][N + 1], bt[N + 1];
int PIV[N + 1];

void Make_Matrix(int n);
void Make_solution(int n);
void Output_Matrix(int n, int line);
void GaussLeft(int n);
void GaussRight(int n);
void Norm(int n);

int main(void)
{
	int i, j, n = N;
	double t, GFLOPS;
	double tic, toc;

	Make_Matrix(N);
	Make_solution(N);
	/* Output_Matrix is used for testing only  */
	/*	Output_Matrix( N, N );	*/

	tic = omp_get_wtime();

	GaussLeft(N);
	GaussRight(N);

	/* Gaussian elimination for the matrix A

	   int PIV[N+1];

	   GaussLeft( N, A, PIV);

	*/
	toc = omp_get_wtime();
	//Norm(N);

	/* Gaussial elimination for the right-hand side vector b

	   GaussRight( N, A, PIV, b);

	*/

	/* if computed-b = original-x = b then your "Gauss" is correct.
	 */

	t = (double)(toc - tic);
	GFLOPS = (2.0 * n) * n * n / t / 1e+9;
	printf("Gaussian Elimination  Size = %3d, Time (sec) = %.2f, GFLOPS = %.2f\n", N, t, GFLOPS);
}

void Make_Matrix(int n)
{
	/* Make a matrix */
	int i, j;
	srand(3);
	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= n; j++)
		{
			A[i][j] = (double)(rand() % 65535) / 65535.0;
			//At[i][j] = A[i][j];
		}
	}
}

void Make_solution(int n)
{
	/* Set solution */
	int i;
	for (i = 1; i <= n; i++)
	{
		b[i] = (double)(rand() % 65535) / 65535.0;
		//bt[i] = b[i];
	}
}

// void Make_solution(int n)
// {
// 	int i, j;
// 	for (i = 1; i <= n; i++)
// 	{
// 		x[i] = (double)(rand() % 65535) / 65535.0;
// 		/*  x[i]= rand()/(double)RAND_MAX;*/
// 	}
// 	/*
// 		printf("put the values of Vector B\n");
// 		double B[FIELD_SIZE];
// 		for (int i = 0; i < line; i++){
// 			scanf("%lf",&B[i]);
// 		}
// 	*/
// 	/* Make a right-hand side vector */
// 	for (i = 1; i <= n; i++)
// 	{
// 		b[i] = 0.0;
// 		for (j = 1; j <= n; j++)
// 		{
// 			b[i] += A[i][j] * x[j];
// 		}
//         bt[i]=b[i];
// 	}
// }

void Output_Matrix(int n, int line)

{
	int i, j;
	/* Confirm the matrix */
	printf("Matrix A:\n");
	for (int i = 1; i <= line; i++)
	{
		for (int j = 1; j <= line; j++)
		{
			// printf("%s%2d%s%7.4lf%s%2d%s%7.4lf","   A[",i+1,",",j+1,"]",A[i][j]);
			printf("%6.2f", A[i][j]);
		}
		putchar('\n');
	}

	/* Confirm the vector */
	printf("Vetctor x and b:\n");
	for (i = 1; i <= n; i++)
	{
		printf(" %5.2f  --  %5.2f\n", x[i], b[i]);
	}
	printf("\n");
}

void GaussLeft(int n)
{
	int i, j, k;
	int pivnum;
	double piv;
	double temp, alpha;
	for (k = 1; k < n; k++)
	{
		piv = fabs(A[k][k]);
		pivnum = k;
		for (i = k + 1; i <= n; i++)
		{
			if (piv < fabs(A[i][k]))
			{
				piv = fabs(A[i][k]);
				pivnum = i;
			}
		}
		if (pivnum != k)
		{
			for (i = 1; i <= n; i++)
			{
				temp = A[k][i];
				A[k][i] = A[pivnum][i];
				A[pivnum][i] = temp;
			}
			temp = b[k];
			b[k] = b[pivnum];
			b[pivnum] = temp;
		}
		if (piv == 0)
		{
			fputs("faild\n", stderr);
			exit(1);
		}
		PIV[k] = pivnum;
		#pragma omp parallel for private(j, alpha)
		for (i = k + 1; i <= n; i++)
		{
			alpha = -A[i][k] / A[k][k];
			A[i][k] = alpha;
			for (j = k + 1; j <= n; j++)
			{
				A[i][j] += alpha * A[k][j];
			}
		}
	}
}

void GaussRight(int n)
{
	int k, j;
	double sum, temp;
	for (k = 2; k <= n; k++)
	{
		for (j = 1; j < k; j++)
		{
			b[k] += A[k][j] * b[j];
		}
	}
	for (k = n; k > 0; k--)
	{
		sum = 0;
		for (j = k + 1; j <= n; j++)
		{
			sum += A[k][j] * b[j];
		}
		b[k] = (b[k] - sum) / A[k][k];
	}
	for (k = n - 1; k > 0; k--)
	{
		if (PIV[k] != k)
		{
			temp = b[k];
			b[k] = b[PIV[k]];
			b[PIV[k]] = temp;
		}
	}
}

// void Norm(int n)
// {
// 	double f[n + 1];
// 	int i, j;
// 	double sum1 = 0;
// 	double sum2 = 0;
// 	for (i = 1; i <= n; i++)
// 	{
// 		f[i] = 0;
// 		for (j = 1; j <= n; j++)
// 		{
// 			f[i] += At[i][j] * b[j];
// 		}
// 	}
// 	for (i = 1; i <= n; i++)
// 	{
// 		f[i] -= bt[i];
// 		sum1 += f[i] * f[i];
// 	}
// 	for (i = 1; i <= n; i++)
// 	{
// 		sum2 += bt[i] * bt[i];
// 	}
// 	printf("ノルム: %lf\n", sqrt(sum1 / sum2));
// }
