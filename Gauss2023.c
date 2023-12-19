#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
/* Kadai_2  2023.12. 8. */
#define N 20
#define BUF_SIZE 256
#define LINE_SIZE 21  // 20×20なので
#define FIELD_SIZE 21 // 20×20なので

double A[N + 1][N + 1], x[N + 1], b[N + 1];
int PIV[N + 1];

void Input_Matrix(int K);
void Make_righthand_side(int n);
void Output_Matrix(int n, int line);
void GaussLeft(int n);
void GaussRight(int n);
void Norm(int n);

int main(void)
{
	int i, j, n = N;
	double t, GFLOPS;
	clock_t tic, toc;

	Input_Matrix(N);
	Make_righthand_side(N);

	/* Output_Matrix is used for testing only  */

	tic = clock();
	GaussLeft(N);
	GaussRight(N);

	/* Gaussian elimination for the matrix A

	   int PIV[N+1];

	   GaussLeft( N, A, PIV);

	*/
	toc = clock();
	//Output_Matrix(N, N);
	Norm(N);

	/* Gaussial elimination for the right-hand side vector b

	   GaussRight( N, A, PIV, b);

	*/

	/* if computed-b = original-x = b then your "Gauss" is correct.
	 */

	t = (double)(toc - tic) / CLOCKS_PER_SEC;
	GFLOPS = (2.0 * n) * n * n / t / 1e+9;
	printf("Gaussian Elimination  Size = %3d, Time (sec) = %.2f, GFLOPS = %.2f\n", N, t, GFLOPS);
}

void Input_Matrix(int K)
{
	/* Only for (LINE_SIZE-1)*(FIELD_SIZE-1) */
	char buf[BUF_SIZE];
	int line;

	FILE *fp = fopen("SparseMatrixB.csv", "r");
	if (fp == NULL)
		return;

	/*  for (line = 0; fgets(buf, BUF_SIZE, fp); line++) ;
		if (line > LINE_SIZE) { return 2; }
		rewind(fp); */

	//
	for (int i = 1; fgets(buf, BUF_SIZE, fp); i++)
	{
		char *p = buf, *sep;
		int j = 0;
		if (line > LINE_SIZE)
		{
			return;
		}
		while (j < FIELD_SIZE)
		{
			A[i][j + 1] = strtod(p, &sep);
			if (++j >= FIELD_SIZE)
			{
				puts("too many fields");
				return;
			}
			if (*sep != ',')
				break;
			p = sep + 1;
		}
	}
	fclose(fp);
}

void Make_righthand_side(int n)
/* Set solution and Right-hand */
{
	int i, j;
	for (i = 1; i <= n; i++)
	{
		x[i] = 1.0;
		/*  x[i]= rand()/(double)RAND_MAX;*/
	}
	/*
		printf("put the values of Vector B\n");
		double B[FIELD_SIZE];
		for (int i = 0; i < line; i++){
			scanf("%lf",&B[i]);
		}
	*/
	/* Make a right-hand side vector */
	for (i = 1; i <= n; i++)
	{
		b[i] = 0.0;
		for (j = 1; j <= n; j++)
		{
			b[i] += A[i][j] * x[j];
		}
	}
}


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
		printf(" %5.2f  --  %lf\n", x[i], b[i]);
	}
	printf("\n");
}

void GaussLeft(int n)
{
	int i, j, k;
	int piv, pivnum;
	double temp, alpha;
	for (k = 1; k < n; k++)
	{
		piv = abs(A[k][k]);
		pivnum = k;
		for (i = k + 1; i <= n; i++)
		{
			if (piv < abs(A[i][k]))
			{
				piv = abs(A[i][k]);
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
		PIV[k] = pivnum;
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

void Norm(int n)
{
	Input_Matrix(n);
	double f[n + 1];
	double c[n + 1];
	int i, j;
	double sum1 = 0, sum2 = 0;
	for (i = 1; i <= n; i++)
	{
		f[i] = 0;
		for (j = 1; j <= n; j++)
		{
			f[i] += A[i][j] * b[j];
		}
	}
	for (i = 1; i <= n; i++)
	{
		c[i] = 0.0;
		for (j = 1; j <= n; j++)
		{
			c[i] += A[i][j] * x[j];
		}
	}
	for (i = 1; i <= n; i++)
	for(i=1;i<=n; i++){
		c[i]=0;
	}
	c[n]=1;
	{
		f[i] -= c[i];
		sum1 += f[i] * f[i];
	}
	for (i = 1; i <= n; i++)
	{
		sum2 += c[i] * c[i];
	}
	printf("ノルム: %lf\n", sqrt(sum1 / sum2));
}
