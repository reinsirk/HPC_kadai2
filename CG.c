/*
1W213026 浦浪 英俊 2024/1/12
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
/* Coojugate Gradient Method  2024. 1.11.  */
#define N 20
#define BUF_SIZE 256
#define LINE_SIZE 21  // 20×20なので
#define FIELD_SIZE 21 // 20×20なので
#define eps 1.0e-13

double A[N + 1][N + 1], x[N + 1], b[N + 1], x_save[N + 1];
double a[(N + 1) * (N + 1)];
int Column_Index[(N + 1) * (N + 1)], Row_Ptr[N + 2];
int PIV[N + 1];

void Input_Matrix(int K);
void Make_righthand_side(int n);
void Output_Matrix(int n, int line);
void ResidualNorm(int n);
void ErrorNorm(int n);
void GaussLeft(int n);
void GaussRight(int n);

void CRS(int n, double *a, int *Column_Index, int *Row_Ptr);
void CRS_multiple(int n, double *y, double *x, double *a, int *Column_Index, int *Row_Ptr);
void A_multiple(int n, double *y, double *x, double *a, int *Column_Index, int *Row_Ptr);
void cg_method(int n, double *x, double *b, double *a, int *Column_Index, int *Row_Ptr);

int main(void)
{
  int i, j, n;
  double t;
  clock_t tic, toc;

  Input_Matrix(N);
  Make_righthand_side(N);
  /* Output_Matrix is used for testing only  */
  Output_Matrix(N, N);

  CRS(N, a, Column_Index, Row_Ptr);

  for (i = 1; i <= N; i++)
  {
    x[i] = 0.0;
  }
  /* tic = clock(); */
  cg_method(N, x, b, a, Column_Index, Row_Ptr);
  GaussLeft(N);
  GaussRight(N);
  ErrorNorm(N);

  /* toc = clock(); */

  /* Please confirm the correctness of computed value. */

  /* t = (double)(toc - tic)/CLOCKS_PER_SEC;
  printf("Conjugate Method  Size = %3d, Time (sec) = %.2f\n", N, t); */
}

void CRS(int n, double *a, int *Column_Index, int *Row_Ptr)
{
  int i, j;
  int sparse_loc = 1;

  for (i = 1; i <= n; i++)
  {
    Row_Ptr[i] = sparse_loc;
    for (j = 1; j <= n; j++)
    {
      if (A[i][j] != 0.0)
      {
        a[sparse_loc] = A[i][j];
        Column_Index[sparse_loc] = j;
        sparse_loc++;
      }
    }

    /* n 行 n 列の密行列から非ゼロ要素を抜き出して CRS 形式に格納する　　*/
  }
  Row_Ptr[i] = sparse_loc;
}

void CRS_multiple(int n, double *y, double *x, double *a, int *Column_Index, int *Row_Ptr)
{
  int i, j, jj;
  double s;

  for (i = 1; i <= n; i++)
  {
    s = 0.0;
    for (j = Row_Ptr[i]; j < Row_Ptr[i + 1]; j++)
    {
      jj = Column_Index[j];
      s += a[j] * x[jj];
    }
    /* この部分を作成し、 A_multiple を CRS_multiple に変更　*/

    y[i] = s;
  }
}

void A_multiple(int n, double *y, double *x, double *a, int *Column_Index, int *Row_Ptr)
{
  int i, j;
  double s;

  for (i = 1; i <= n; i++)
  {
    s = 0.0;
    for (j = 1; j <= n; j++)
    {
      s = s + A[i][j] * x[j];
    }
    y[i] = s;
  }
}

void cg_method(int n, double *x, double *b, double *a, int *Column_Index, int *Row_Ptr)
{
  double p[N + 1], r[N + 1], Ax[N + 1], Ap[N + 1];
  double alpha, beta, bn, rr, rr_new, apr;
  int i, count = 0;

  CRS_multiple(n, Ax, x, a, Column_Index, Row_Ptr);

  bn = 0.0;
  rr = 0.0;
  for (i = 1; i <= n; i++)
  {
    p[i] = b[i] - Ax[i];
    r[i] = p[i];
    rr = rr + r[i] * r[i];
    bn = bn + b[i] * b[i];
  }

  while (1)
  {
    count++;

    CRS_multiple(n, Ap, p, a, Column_Index, Row_Ptr);
    apr = 0.0; /* inner_product(n, p, Ap); */
    for (i = 1; i <= n; i++)
    {
      apr = apr + p[i] * Ap[i];
    }
    alpha = rr / apr;

    rr_new = 0.0;
    for (i = 1; i <= n; i++)
    {
      x[i] += alpha * p[i];
      r[i] -= alpha * Ap[i];
      rr_new = rr_new + r[i] * r[i];
    }

    printf("Iteration    %3d: %2.3e\n", count, sqrt(rr_new / bn));

    if (eps > sqrt(rr_new / bn))
      break;

    beta = rr_new / rr;
    for (i = 1; i <= n; i++)
    {
      p[i] = r[i] + beta * p[i];
    }
    rr = rr_new;
  }

  CRS_multiple(n, Ax, x, a, Column_Index, Row_Ptr);
  rr = 0.0;
  for (i = 1; i <= n; i++)
  {
    r[i] = b[i] - Ax[i];
    rr = rr + r[i] * r[i]; /* sqrt(inner_product(n, r, r) */
  }

  printf("Final Residual :  %2.3e\n", sqrt(rr / bn));
  printf("\n");

  printf("Vetctor x and xsave:\n");
  for (i = 1; i <= n; i++)
  {
    printf(" %5.15f  --  %5.15f\n", x[i], x_save[i]);
  }
  printf("\n");
  ResidualNorm(n);
}

void Input_Matrix(int K)
{
  /* Only for (LINE_SIZE-1)*(FIELD_SIZE-1) */
  char buf[BUF_SIZE];
  int line;

  FILE *fp = fopen("SparseMatrixA.csv", "r");
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
    x_save[i] = x[i];
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
    printf(" %5.2f  --  %5.2f\n", x[i], b[i]);
  }
  printf("\n");
}

void ResidualNorm(int n)
{
  double f, norm;
  int i, j;
  double sum1 = 0.0;
  double sum2 = 0.0;
  for (i = 1; i <= n; i++)
  {
    f = 0;
    for (j = 1; j <= n; j++)
    {
      f += A[i][j] * x[j];
    }
    f -= b[i];
    sum1 += pow(f, 2);
    // sum2 += pow(b[i], 2);
  }
  // norm = sqrt(sum1 / sum2);
  norm = sqrt(sum1);
  printf("Residual Norm(残差ノルム): %2.3e\n", norm);
}

void ErrorNorm(int n)
{
  double f, norm;
  int i, j;
  double sum1 = 0.0;
  double sum2 = 0.0;
  for (i = 1; i <= n; i++)
  {
    f = x[i] - b[i];
    sum1 += pow(f, 2);
    // sum2 += pow(b[i], 2);
  }
  // norm = sqrt(sum1 / sum2);
  norm = sqrt(sum1);
  printf("Error Norm(CG法の解とガウスの消去法の解の誤差ノルム): %2.3e\n", norm);
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
      temp = fabs(A[i][k]);
      if (piv < temp)
      {
        piv = temp;
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
    }
    if (piv == 0.0)
    {
      fputs("faild\n", stderr);
      exit(1);
    }
    PIV[k] = pivnum;
    for (i = k + 1; i <= n; i++)
    {
      alpha = -1 * (A[i][k] / A[k][k]);
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
  for (k = 1; k < n; k++)
  {
    if (PIV[k] != k)
    {
      temp = b[k];
      b[k] = b[PIV[k]];
      b[PIV[k]] = temp;
    }
  }
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
}