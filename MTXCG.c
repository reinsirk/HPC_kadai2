#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
/* Coojugate Gradient Method  2024. 1.17.  */
#define N 4000
#define NN (N + 1) * (N + 1) /* Enough. But better to allocate dynamicaly using nnz */
#define eps 1.0e-5

double x[N + 1], b[N + 1], x_save[N + 1], a[NN];
int Column_Index[(N + 1) * (N + 1)], Row_Ptr[N + 2];

void MTXCRS(int *n, int *nnz, double *a, int *Column_Index, int *Row_Ptr, int *ecode);
void CRS_multiple(int n, double *y, double *x, double *a, int *Column_Index, int *Row_Ptr);
void CRS_multiple_sym(int n, double *y, double *x, double *a, int *Column_Index, int *Row_Ptr);
void cg_method(int n, double *x, double *b, double *a, int *Column_Index, int *Row_Ptr, int *nct, double *res_fin);
void Make_righthand_side(int n);

int main(void)
{
  int n, nnz, nct, i, ecode;
  double t, res_fin;
  clock_t tic, toc;

  /* Set Lower part of Symmetric matrix as CRS format and solution, right-hand side */
  MTXCRS(&n, &nnz, a, Column_Index, Row_Ptr, &ecode);
  Make_righthand_side(n);

  while (0)
  { /* Set 1 for Debug */
    printf("n: %d   nnz : %d\n", n, nnz);
    printf("Content of Row_Ptr:\n");
    for (i = 1; i <= n + 1; i++)
    {
      printf("  %d  -  %d\n", i, Row_Ptr[i]);
    }
    printf("Cooum_Index and Value:\n");
    for (i = 1; i <= Row_Ptr[n + 1] - 1; i++)
    {
      printf(" %d Col: %d - Val: %f\n", i, Column_Index[i], a[i]);
    }
    printf("Vetctor x and b:\n");
    for (i = 1; i <= n; i++)
    {
      printf(" %5.2f  --  %5.2f\n", x[i], b[i]);
    }
    printf("\n");
    break;
  }
  /* Set an initial solution. Do not use the real solution. */
  for (i = 1; i <= N; i++)
  {
    x[i] = 0.0;
  }

  tic = clock();
  /* CG method including confirmation of computed solution. */
  cg_method(n, x, b, a, Column_Index, Row_Ptr, &nct, &res_fin);
  toc = clock();
  t = (double)(toc - tic) / CLOCKS_PER_SEC;

  printf("# of Iteration     =  %8d\n", nct);
  printf("Relative Residual  =  %2.3e\n", res_fin);
  printf("Computaion time    = %.2f\n", t);
  printf("\n");
}
/* 　対称行列の上／下三角部分が与えられた場合　　*/
void CRS_multiple_sym(int n, double *y, double *x, double *a, int *Column_Index, int *Row_Ptr)
{
  int i, j, jj;
  double s;

  for (i = 1; i <= n; i++)
  {
    y[i] = 0.0;
  }

  for (i = 1; i <= n; i++)
  {
    s = 0.0;
    for (j = Row_Ptr[i]; j < Row_Ptr[i + 1]; j++)
    {
      jj = Column_Index[j];
      if (i == jj)
      {
        s += a[j] * x[jj]; /* 対角:： A[i][i] に対して */
      }
      else
      {
        s += a[j] * x[jj];    /* 非対角： y[i] += A[i][j]*x[j] に対する式　*/
        y[jj] += a[j] * x[i]; /* 対称な場所　 　y[j] += A[j][i]*x[i] に対する式　*/
      }
    }
    y[i] += s;
  }
}
/* 　一般的なＣＲＳ形式での行列ベクトル積（前回の解答例）　*/
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
      /*    if ( jj <= i ) {   */ /* 下三角部分 = if 文の間だけを有効にして完全な行列ベクトル積を構築するのと同値　　　　*/
      s += a[j] * x[jj];
      /*    }    */
    }
    y[i] = s;
  }
}

void cg_method(int n, double *x, double *b, double *a, int *Column_Index, int *Row_Ptr, int *nct, double *res_fin)
{
  double p[N + 1], r[N + 1], Ax[N + 1], Ap[N + 1];
  double alpha, beta, bn, rr, rr_new, apr;
  int i, count = 0;

  CRS_multiple_sym(n, Ax, x, a, Column_Index, Row_Ptr);

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

    CRS_multiple_sym(n, Ap, p, a, Column_Index, Row_Ptr);
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
    {
      break;
    }
    beta = rr_new / rr;
    for (i = 1; i <= n; i++)
    {
      p[i] = r[i] + beta * p[i];
    }
    rr = rr_new;
  }

  CRS_multiple_sym(n, Ax, x, a, Column_Index, Row_Ptr);
  rr = 0.0;
  for (i = 1; i <= n; i++)
  {
    r[i] = b[i] - Ax[i];
    rr = rr + r[i] * r[i]; /* sqrt(inner_product(n, r, r) */
  }
  *nct = count;
  *res_fin = sqrt(rr / bn);

  printf("Final Residual =  %2.3e\n", *res_fin);
  printf("\n");
  printf("Vetctor x and xsave:\n");
  for (i = 1; i <= n; i++)
  {
    printf(" %5.15f  --  %5.15f  :  %2.3e\n", x[i], x_save[i], fabs(x[i] - x_save[i]));
  }
  printf("\n");
}

/* Convert MatrixMarket format with row-wise order to CRS format  */
void MTXCRS(int *n, int *nnz, double *a, int *Column_Index, int *Row_Ptr, int *ecode)
{
  int i, j, na, ma, nnza, loop, nn, k;
  double value;
  char line[100];
  long int position;
  int sparse_loc = 1;
  int ilast = 0;
  na=2;
  // ファイルを開く
  FILE *fp;
  fp = fopen("bcsstk15.mtx", "r");
  if (fp == NULL)
  {
    perror("ファイルを開けませんでした");
    *ecode = 1;
    exit(1);
  }
  // 先頭が '%'の行を読み飛ばす
  while (fgets(line, sizeof(line), fp))
  {
    /* line[100] = '\0';
    printf("Input: %s\n",line); */
    if (line[0] != '%')
    {
      break;
    }
    position = ftell(fp);
  }
  // 3つの整数を読み取る
  fseek(fp, position, SEEK_SET);
  if (fscanf(fp, "%d %d %d\n", &na, &ma, &nnza) != 3)
  {
    perror("n, m, nnz を読み取れませんでした");
    *ecode = 1;
    exit(1);
  }
  /* printf("%d %d %d\n", nn, m, nnz); */
  *n = na;
  *nnz = nnza;
  loop = nnza;
  position = ftell(fp);
  // 行、列、値の組み合わせを読み取る
  // データ区切りに , がないときは要修正
  for (k=1;k<=na;k++){
  Row_Ptr[k] = sparse_loc;
  fseek(fp, position, SEEK_SET);
  loop = nnza;
  while (loop)
  {
    if (fscanf(fp, "%d %d %lf\n", &i, &j, &value) != 3)
    {
      perror("データを読み取れませんでした");
      *ecode = 1;
      exit(1);
    }
    if(k<j){
      break;
    }
    /* CRS : input must be a row-wise ordered */
    /*if (i != ilast)
    {
      Row_Ptr[i] = sparse_loc;
    }
    Column_Index[sparse_loc] = j;
    a[sparse_loc] = value;
    sparse_loc = sparse_loc + 1;
    ilast = i;*/
    /*        printf("%d, %d, %f\n", i, j, value); */
    if(i==k){
    Column_Index[sparse_loc] = j;
    a[sparse_loc] = value;
    sparse_loc = sparse_loc + 1;
    }
    loop = loop - 1;
  }
  //printf("%d\n",k);
  }
  Row_Ptr[k] = sparse_loc;  
  fclose(fp);
  *ecode = 0;
  return;
}

void Make_righthand_side(int n)
/* Set solution and Right-hand */
{
  int i;
  for (i = 1; i <= n; i++)
  {
    x[i] = 1.0;
    x_save[i] = x[i];
    /*  x[i]= rand()/(double)RAND_MAX;*/
  }
  /* Make a right-hand side vector */
  CRS_multiple_sym(n, b, x, a, Column_Index, Row_Ptr);
}

