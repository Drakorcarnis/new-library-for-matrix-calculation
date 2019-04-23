#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <complex.h>
#include <immintrin.h>
#include <omp.h>
#include "matrix.h"
#include "matrix_tools.h"

#define MAX_RANK 10000
#define NUM_POOLS 2

#define ROW_MODE    0
#define COLUMN_MODE 1

typedef struct {
    int nb_perm;
    matrix_t *P;
    matrix_t *L;
    matrix_t *U;
} plu_t;

typedef struct {
    double      cached_matrix[MAX_RANK*MAX_RANK];
    void *      cached_matrix_ptr;
    int         cached_matrix_mode;
    int         index;
}pool_t;

pool_t pools[NUM_POOLS]; 
int pool_index = 0; 

static int sanity_check(const void *pointer, const char *function_name);
static int square_check(const matrix_t *matrix, const char *function_name);
static int symetry_check(const matrix_t *matrix, const char *function_name);

static plu_t * plu_create(unsigned int rank);
static void plu_free(plu_t *plu);
static plu_t * matrix_plu_f(const matrix_t *matrix);
static double complex * matrix_cholesky_f(const matrix_t *matrix);

static int sanity_check(const void *pointer, const char *function_name)
{
    if(!pointer){
        fprintf(stderr, "%s: NULL pointer\n",function_name);
        return 0;
    }
    return 1;
}

static int square_check(const matrix_t *matrix, const char *function_name)
{
    if((matrix->rows != matrix->columns)){
        fprintf(stderr, "%s: not square matrix\n",function_name);
        return 0;
    }
    return 1;
}

static int symetry_check(const matrix_t *matrix, const char *function_name)
{
    for (unsigned int i = 0; i < matrix->rows; i++) {
        for (unsigned int j = 0; j < i; j++) {
            if(matrix->coeff[i][j] != matrix->coeff[j][i]){
                fprintf(stderr, "%s: not symetric matrix\n",function_name);
                return 0;
            }
        }
    }
    return 1;
}

static double * stackify_matrix(const matrix_t *matrix, int mode)
{
    double * ret = NULL;
    for (unsigned int i = 0; i < NUM_POOLS; i++) {
        if (pools[i].cached_matrix_ptr == matrix && pools[i].cached_matrix_mode == mode){
            return pools[i].cached_matrix;
        }
    }
    for (unsigned int i = 0; i < NUM_POOLS; i++) {
        if (pools[i].cached_matrix_ptr == NULL){
            ret = pools[i].cached_matrix;
            pools[i].cached_matrix_ptr = (void *) matrix;
            pools[i].cached_matrix_mode = mode;
            break;
        }
    }
    if (ret == NULL){
        pools[pool_index].cached_matrix_mode = mode;
        pools[pool_index].cached_matrix_ptr = (void *) matrix;
        ret = pools[pool_index].cached_matrix;
        pool_index = (pool_index + 1) % NUM_POOLS;
    }
    unsigned int n = matrix->rows;
    unsigned int m = matrix->columns;
    if (n>MAX_RANK){
        fprintf(stderr,"Rank %u > MAX_RANK\n", n);
        return NULL;
    }
    if (m>MAX_RANK){
        fprintf(stderr,"Rank %u > MAX_RANK\n", m);
        return NULL;
    }
    if(mode == ROW_MODE){
        for (unsigned int i = 0; i < n; i++) {
            for (unsigned int j = 0; j < m; j++) {
                ret[i * n + j] = matrix->coeff[i][j];
            }
        }
    } else if(mode == COLUMN_MODE){
        for (unsigned int i = 0; i < n; i++) {
            for (unsigned int j = 0; j < m; j++) {
                ret[j * m + i] = matrix->coeff[i][j];
            }
        }
    } else {
        fprintf(stderr,"Invalid mode %d\n", mode);
        return NULL;
    }
    return ret;
}

static plu_t * plu_create(unsigned int rank){
    plu_t *plu = malloc(sizeof(plu_t));
    plu->P = matrix_identity(rank);
    plu->L = matrix_create(rank, rank);
    plu->U = matrix_identity(rank);
    plu->nb_perm = 0;
    return plu;
}

static void plu_free(plu_t *plu)
{
    if(!sanity_check(plu, __func__))return; 
    matrix_free(plu->P);
    matrix_free(plu->L);
    matrix_free(plu->U);
    free(plu);   
}
// Matrix creation functions

matrix_t * matrix_create(unsigned int rows, unsigned int columns)
{
    unsigned int i = 0;
    matrix_t *matrix = malloc(sizeof(matrix_t));
    if (!matrix) goto failed_matrix;
    matrix->rows = rows;
    matrix->columns = columns;
    matrix->coeff = malloc(rows*sizeof(double *));
    if (!matrix->coeff) goto failed_coeff;
    for (; i < rows; i++){
        matrix->coeff[i] = calloc(columns, sizeof(double));
        if (!matrix->coeff[i]) goto failed_coeff_elt;
    }
    return matrix;
failed_coeff_elt:
    for (unsigned int j = 0; j < i; j++)free(matrix->coeff[j]);
    free(matrix->coeff);
failed_coeff:
    free(matrix);
failed_matrix:
    perror("alloc");
    return NULL;
}

matrix_t * matrix_identity(unsigned int n)
{
    matrix_t * matrix = matrix_create(n, n);
    for (unsigned int i = 0; i < n; i++) matrix->coeff[i][i] = 1;
    return matrix;
}

matrix_t * matrix_permutation(unsigned int line1, unsigned int line2, unsigned int n)
{
    double tmp;
    matrix_t * matrix = matrix_identity(n);
    for (unsigned int i = 0; i < n; i++){
        tmp = matrix->coeff[line1][i];
        matrix->coeff[line1][i] = matrix->coeff[line2][i];
        matrix->coeff[line2][i] = tmp;
    }
    matrix_t * transpose = matrix_transp_f(matrix);
    matrix_free(matrix);
    return(transpose);
}

matrix_t * matrix_copy(const matrix_t *matrix)
{
    matrix_t *copy = matrix_create(matrix->rows, matrix->columns);
    for (unsigned int i = 0; i < matrix->rows; i++) {
        for (unsigned int j = 0; j < matrix->columns; j++) {
            copy->coeff[i][j] = matrix->coeff[i][j];
        }
    }
    return copy;
}

void matrix_free(matrix_t *matrix)
{
    if(!sanity_check(matrix, __func__))return; 
    for (unsigned int i = 0; i < matrix->rows; i++) {
        free(matrix->coeff[i]);
    }
    free(matrix->coeff);
    free(matrix); 
    matrix = NULL;
}

// Matrix computation functions

// Basic operations
matrix_t * matrix_transp_f(const matrix_t *matrix)
{
    if(!sanity_check((void *)matrix, __func__))return NULL; 
    matrix_t *transpose_matrix = matrix_create(matrix->columns, matrix->rows);
    for (unsigned int i = 0; i < transpose_matrix->rows; i++) {
        for (unsigned int j = 0; j < transpose_matrix->columns; j++)transpose_matrix->coeff[i][j]=matrix->coeff[j][i];
    }
    return transpose_matrix;
}

matrix_t * matrix_add_f(const matrix_t *matrix1, const matrix_t *matrix2)
{  
    if(!sanity_check((void *)matrix1, __func__))return NULL; 
    if(!sanity_check((void *)matrix2, __func__))return NULL; 
    if((matrix1->rows != matrix2->rows) || (matrix1->columns != matrix2->columns)){
        fprintf(stderr, "%s: not addable matrix ((matrix1->rows != matrix2->rows) || (matrix1->columns != matrix2->columns))\n", __func__);
        return NULL;
    }
    matrix_t *add_matrix = matrix_create(matrix1->rows, matrix1->columns);
    for (unsigned int i = 0; i < add_matrix->rows; i++) {
        for (unsigned int j = 0; j < add_matrix->columns; j++) {
            add_matrix->coeff[i][j] = matrix1->coeff[i][j] + matrix2->coeff[i][j];
        }
    }
    return add_matrix;
}

matrix_t * matrix_mult_scalar_f(const matrix_t *matrix, double lambda)
{
    if(!sanity_check((void *)matrix, __func__))return NULL;
    matrix_t *mult_matrix = matrix_create(matrix->rows, matrix->columns);
    for (unsigned int i = 0; i < matrix->rows; i++) {
        for (unsigned int j = 0; j < matrix->columns; j++) {
            mult_matrix->coeff[i][j]=lambda * matrix->coeff[i][j];
        }
    }
    return mult_matrix;
}

static double avx_mult(int n, double *x, double *y)
{
	int i;
    int n8 = n>>3<<3;
	__m256d vs1, vs2;
	double s, t[4];
	vs1 = _mm256_setzero_pd();
	vs2 = _mm256_setzero_pd();
	for (i = 0; i < n8; i += 8) {
		__m256d vx1, vx2, vy1, vy2;
		vx1 = _mm256_loadu_pd(&x[i]);
		vx2 = _mm256_loadu_pd(&x[i+4]);
		vy1 = _mm256_loadu_pd(&y[i]);
		vy2 = _mm256_loadu_pd(&y[i+4]);
		vs1 = _mm256_add_pd(vs1, _mm256_mul_pd(vx1, vy1));
		vs2 = _mm256_add_pd(vs2, _mm256_mul_pd(vx2, vy2));
	}
	for (s = 0.0f; i < n; ++i){
        s += x[i] * y[i];
    }
    
	_mm256_storeu_pd(t, vs1);
	s += t[0] + t[1] + t[2] + t[3];
	_mm256_storeu_pd(t, vs2);
	s += t[0] + t[1] + t[2] + t[3];
	return s;
}

matrix_t * matrix_mult_f(const matrix_t *matrix1, const matrix_t *matrix2)
{
    if(!sanity_check((void *)matrix1, __func__))return NULL; 
    if(!sanity_check((void *)matrix2, __func__))return NULL; 
    if((matrix2->rows != matrix1->columns)){
        fprintf(stderr, "%s: not multiplicable matrix (matrix2->rows != matrix1->columns)\n", __func__);
        return NULL;
    }
   	unsigned int i, j, n = matrix1->rows, m = matrix2->columns;
    matrix_t *mult = matrix_create(matrix1->columns, matrix2->columns);
    double *rows = stackify_matrix(matrix1, ROW_MODE);
    if(!sanity_check((void *)rows, __func__))return NULL; 
    double *columns = stackify_matrix(matrix2, COLUMN_MODE);
    if(!sanity_check((void *)columns, __func__))return NULL; 
    double ** coeff = mult->coeff;
    #pragma omp parallel shared(coeff) private(i, j) num_threads(8)
	{
        #pragma omp for schedule(static)
        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                coeff[i][j] = avx_mult(n, &rows[i*n], &columns[j*m]);
            }
        }
    }
	return mult;
}

matrix_t * matrix_pow_f(const matrix_t *matrix, int pow)
{
    if(!sanity_check((void *)matrix, __func__))return NULL;  
    if(!square_check(matrix, __func__))return NULL; 
    matrix_t *pow_matrix = matrix_copy(matrix);
    matrix_t *tmp_matrix;
    for (int i = 0; i < pow-1; i++) {
        tmp_matrix = matrix_mult_f(pow_matrix, matrix);
        for (unsigned int j = 0; j < pow_matrix->rows; j++) {
            for (unsigned int k = 0; k < pow_matrix->columns; k++) {
                pow_matrix->coeff[j][k] = tmp_matrix->coeff[j][k];
            }
        }
        matrix_free(tmp_matrix);
    }
    return pow_matrix;
}

matrix_t * matrix_shrink_f(const matrix_t *matrix, unsigned int skipped_row, unsigned int skipped_column)
{
    if(!sanity_check((void *)matrix, __func__))return NULL;
    unsigned int l = 0;
    matrix_t *shrink_matrix = matrix_create(matrix->rows-1, matrix->columns-1);
    for (unsigned int i = 0; i < shrink_matrix->rows; i++) {
        if (l == skipped_row) l++;
        unsigned int k=0;
        for (unsigned int j = 0; j < shrink_matrix->columns; j++) {
            if (k == skipped_column) k++;
            shrink_matrix->coeff[i][j]=matrix->coeff[l][k++];
        }
        l++;
    }
    return shrink_matrix;
}

matrix_t * matrix_solve_diag_inf(const matrix_t *A, const matrix_t *B)
{
    if(!sanity_check((void *)A, __func__))return NULL;
    if(!square_check(A, __func__))return NULL; 
    unsigned int n = A->rows;
    unsigned int m = B->columns;
    matrix_t *X = matrix_create(n, n);
    for (unsigned int i = 0; i < m; i++)
    {
        for (unsigned int j = 0; j < n; j++)
        {
            X->coeff[j][i] = B->coeff[j][i];
            for (unsigned int k = 0; k < j; k++)X->coeff[j][i] =  X->coeff[j][i] - X->coeff[k][i] * A->coeff[j][k];
            X->coeff[j][i] = X->coeff[j][i] / A->coeff[j][j];
        }
    }
    return(X);
}

matrix_t * matrix_solve_diag_sup(const matrix_t *A, const matrix_t *B)
{
    if(!sanity_check((void *)A, __func__))return NULL;
    if(!square_check(A, __func__))return NULL; 
    unsigned int n = A->rows;
    unsigned int m = B->columns;
    matrix_t *X = matrix_create(n, n);
    for (unsigned int i = 0; i < m; i++){
        for (int j = n - 1; j >= 0; j--){
            X->coeff[j][i] = B->coeff[j][i];
            for (int k = n - 1; k > j; k--)X->coeff[j][i] =  X->coeff[j][i] - X->coeff[k][i] * A->coeff[j][k];
            X->coeff[j][i] = X->coeff[j][i] / A->coeff[j][j];
        }
    }
    return(X);
}
 
// Methods based upon raw determinant calculation. For fun only. Do never use them, cuz you've NO reason to use them. Really.
double matrix_det_raw_f(const matrix_t *matrix)
{
    double det = 0;
    double sign = 0;
    matrix_t *shrinked_matrix = NULL;
    if(!sanity_check((void *)matrix, __func__))return 0;
    if(matrix->columns != matrix->rows) return 0;
    if(matrix->columns == 1)return matrix->coeff[0][0];
    for (unsigned int i = 0; i < matrix->columns; i++) 
    {
        shrinked_matrix=matrix_shrink_f(matrix, 0, i);
        sign = (int)pow(-1.0,(double)i);
        det += sign * matrix->coeff[0][i] * matrix_det_raw_f(shrinked_matrix);
        matrix_free(shrinked_matrix);
    }
    return det;
}

matrix_t * matrix_inverse_raw_f(const matrix_t *matrix)
{
    if(!sanity_check((void *)matrix, __func__))return NULL;
    double det = matrix_det_raw_f(matrix);
    double one = 1;
    if(!det){
        fprintf(stderr, "%s: not inversible matrix (|M| = 0)\n", __func__);
        return NULL;
    }
    matrix_t *compl_matrix = matrix_comp_f(matrix);
    matrix_t *inverse_matrix = matrix_mult_scalar_f(compl_matrix, one / det);
    matrix_free(compl_matrix);
    return(inverse_matrix);
}

matrix_t * matrix_solve_raw_f(const matrix_t *A, const matrix_t *B)
{
    if(!sanity_check((void *)A, __func__))return NULL;
    if(!sanity_check((void *)B, __func__))return NULL;
    if(!square_check(A, __func__))return NULL;
    matrix_t *inverse_matrix = matrix_inverse_raw_f(A);
    if(!sanity_check((void *)inverse_matrix, __func__))return NULL;
    matrix_t *X = matrix_mult_f(inverse_matrix, B);
    matrix_free(inverse_matrix);
    return(X);
}

matrix_t * matrix_com_f(const matrix_t *matrix)
{
    if(!sanity_check((void *)matrix, __func__))return NULL;
    if(!square_check(matrix, __func__))return NULL; 
    double sign = 0;
    matrix_t *shrinked_matrix = NULL;
    matrix_t *co_matrix = matrix_create(matrix->rows, matrix->columns);
    for (unsigned int i = 0; i < co_matrix->rows; i++){
        for (unsigned int j = 0; j < co_matrix->columns; j++){
            shrinked_matrix=matrix_shrink_f(matrix, i, j);
            sign = (int)pow(-1.0,(double)(i+j));
            co_matrix->coeff[i][j] = sign * matrix_det_raw_f(shrinked_matrix);
            matrix_free(shrinked_matrix);
        }
    }
    return co_matrix;
}

matrix_t * matrix_comp_f(const matrix_t *matrix)
{
    if(!sanity_check((void *)matrix, __func__))return NULL;
    if(!square_check(matrix, __func__))return NULL; 
    matrix_t *co_matrix = matrix_com_f(matrix);
    matrix_t *compl_matrix = matrix_transp_f(co_matrix);
    matrix_free(co_matrix);
    return compl_matrix;
}

// Highly optimized fast methods based upon PLU decomposition.
static plu_t * matrix_plu_f(const matrix_t *matrix)
{
    if(!sanity_check((void *)matrix, __func__))return NULL;
    if(!square_check(matrix, __func__))return NULL;
    unsigned int p, n = matrix->rows;
    plu_t *plu = plu_create(n); 
    matrix_t *A = matrix_copy(matrix), *L = plu->L, *U = plu->U;
    matrix_t *perm, *perm_mult;
    double sum, tmp;
    for (unsigned int i = 0; i < n; i++)
    {
        for (unsigned int j = i; j < n; j++){
            sum = 0;
            for (unsigned int k = 0; k < i; k++)sum += L->coeff[j][k] * U->coeff[k][i];
            L->coeff[j][i] = A->coeff[j][i] - sum;
        }
        // permutations de lignes si pivot nul 
        p = i;
        while((p < n-1) && L->coeff[p][i] == 0 && ++p);
        if (p != i)
        {
            perm = matrix_permutation(p, i, n);
            perm_mult = matrix_mult_f(perm, plu->P);
            matrix_free(perm);
            matrix_free(plu->P);
            plu->P = perm_mult;
            plu->nb_perm+=1;
            for (unsigned int j = 0; j < n; j++){
                tmp = A->coeff[i][j];
                A->coeff[i][j] = A->coeff[p][j];
                A->coeff[p][j] = tmp;
            }
            for (unsigned int j = 0; j < n; j++){
                tmp = L->coeff[i][j];
                L->coeff[i][j] = L->coeff[p][j];
                L->coeff[p][j] = tmp;
            }
        } // Fin des permutations
        for (unsigned int j = i+1; j < n; j++)
        {
            sum = 0;
            for (unsigned int k = 0; k < i; k++){
                sum += L->coeff[i][k] * U->coeff[k][j];
            }  
            U->coeff[i][j] = (A->coeff[i][j] - sum) / L->coeff[i][i];
        }
    }
    matrix_free(A);
    return(plu);  
}

double matrix_det_plu_f(const matrix_t *matrix)
{
    if(!sanity_check((void *)matrix, __func__))return 0;
    if(!square_check(matrix, __func__))return 0; 
    plu_t *plu = matrix_plu_f(matrix);
    double det = pow(-1.0, plu->nb_perm);
    for (unsigned int i=0; i < plu->L->rows; i++)det *= plu->L->coeff[i][i];
    plu_free(plu);
    return det;
}
  
matrix_t * matrix_solve_plu_f(const matrix_t *A, const matrix_t *B)
{
    if(!sanity_check((void *)A, __func__))return NULL;
    if(!square_check(A, __func__))return NULL; 
    plu_t *plu = matrix_plu_f(A);
    matrix_t *transp_perm = matrix_transp_f(plu->P);
    matrix_t *new_B = matrix_mult_f(transp_perm, B);
    matrix_free(transp_perm);
    matrix_t *Z = matrix_solve_diag_inf(plu->L, new_B);
    matrix_free(new_B);
    matrix_t *X = matrix_solve_diag_sup(plu->U, Z);
    matrix_free(Z);
    plu_free(plu);
    matrix_t *ret = matrix_create(A->rows,B->columns);
    for (unsigned int i=0; i < A->rows; i++){
        for (unsigned int j=0; j < B->columns; j++){
            ret->coeff[i][j] = X->coeff[i][j];
        }
    }
    matrix_free(X);
    return(ret);
} 

matrix_t * matrix_inverse_plu_f(const matrix_t *matrix)
{
    if(!sanity_check((void *)matrix, __func__))return NULL;
    if(!square_check(matrix, __func__))return NULL; 
    matrix_t *Id = matrix_identity(matrix->rows);
    matrix_t *matrix_inverse = matrix_solve_plu_f(matrix, Id);
    matrix_free(Id);
    return(matrix_inverse);
} 

static double complex * matrix_cholesky_f(const matrix_t *matrix)
{
    if(!sanity_check((void *)matrix, __func__))return NULL;
    if(!square_check(matrix, __func__))return NULL;
    if(!symetry_check(matrix, __func__))return NULL;
    int n = matrix->rows;
    double complex sum;
    double complex *L = calloc(n*n,sizeof(double complex)); 
    if(!L){
        perror("alloc");
        return NULL;
    }
    for (int i = 0; i < n; i++){
        sum = 0;
        int index = n*i;
        for (int k = 0; k < i; k++)sum += L[index+k] * L[index+k];
        L[n*i+i] = csqrt(matrix->coeff[i][i] - sum);
        for (int j = i+1; j < n; j++){
            sum = 0;
            int index2 = n*j;
            for (int k = 0; k < i; k++)sum += L[index+k] * L[index2+k];
            L[index2+i] = (matrix->coeff[i][j] - sum)/L[index+i];
        }
    }
    return(L);  
} 

double complex ** matrix_solve_diag_inf_comp(double complex *A[], double complex *B[], int n, int m)
{
    double complex **X = calloc(n,sizeof(double complex *));
    if(!X){
        perror("alloc");
        return NULL;
    }
    for (int i = 0; i < n; i++){
        X[i] = calloc(m,sizeof(double complex));
    }
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            X[j][i] = B[j][i];
            for (int k = 0; k < j; k++)X[j][i] =  X[j][i] - X[k][i] * A[j][k];
            X[j][i] = X[j][i] / A[j][j];
        }
    }
    return(X);
}

double complex ** matrix_solve_diag_sup_comp(double complex *A[], double complex *B[], int n, int m)
{
    double complex **X = calloc(n,sizeof(double complex *));
    if(!X){
        perror("alloc");
        return NULL;
    }
    for (int i = 0; i < n; i++){
        X[i] = calloc(m,sizeof(double complex));
    }
    for (int i = 0; i < m; i++){
        for (int j = n - 1; j >= 0; j--){
            X[j][i] = B[j][i];
            for (int k = n - 1; k > j; k--)X[j][i] =  X[j][i] - X[k][i] * A[j][k];
            X[j][i] = X[j][i] / A[j][j];
        }
    }
    return(X);
} 

matrix_t * matrix_solve_cholesky_f(const matrix_t *A, const matrix_t *B)
{
    if(!sanity_check((void *)A, __func__))return NULL;
    if(!square_check(A, __func__))return NULL;
    if(!symetry_check(A, __func__))return NULL;
    int n = A->rows, m = B->columns;
    double complex **comp_B = calloc(n,sizeof(double complex *));
    if(!comp_B){
        perror("alloc");
        return NULL;
    }
    for (int i = 0; i < n; i++){
        comp_B[i] = calloc(m,sizeof(double complex));
    }
    for (int i=0; i < n; i++){
        for (int j=0; j < m; j++){
            comp_B[i][j] = (double complex)B->coeff[i][j];
        }
    }
    double complex *L = matrix_cholesky_f(A);
    double complex **stackL = calloc(n,sizeof(double complex *));
    if(!stackL){
        perror("alloc");
        return NULL;
    }
    double complex **stackLT = calloc(n,sizeof(double complex *));
    if(!stackLT){
        perror("alloc");
        return NULL;
    }
    for (int i = 0; i < n; i++){
        stackL[i] = calloc(n,sizeof(double complex));
        stackLT[i] = calloc(n,sizeof(double complex));
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            stackL[i][j] = L[i*n+j];
            stackLT[j][i] = L[i*n+j];
        }
    }
    free(L);
    double complex **Z = matrix_solve_diag_inf_comp(stackL, comp_B, n, m);
    for (int i = 0; i < n; i++){
        free(comp_B[i]);
        free(stackL[i]);
    }
    free(comp_B);
    free(stackL);
    double complex **X = matrix_solve_diag_sup_comp(stackLT, Z, n, m);
    for (int i = 0; i < n; i++){
        free(Z[i]);
        free(stackLT[i]);
    }
    free(Z);
    free(stackLT);
    matrix_t *ret = matrix_create(n,m);
    for (int i=0; i < n; i++){
        for (int j=0; j < m; j++){
            ret->coeff[i][j] = creal(X[i][j]);
        }
    }
    for (int i = 0; i < n; i++){
        free(X[i]);
    }
    free(X);
    return(ret);
}

matrix_t * matrix_inverse_cholesky_f(const matrix_t *matrix)
{
    if(!sanity_check((void *)matrix, __func__))return NULL;
    if(!square_check(matrix, __func__))return NULL;
    if(!symetry_check(matrix, __func__))return NULL;
    matrix_t *Id = matrix_identity(matrix->rows);
    // matrix_t *matrix_inverse2 = matrix_inverse_plu_f(matrix);

    // matrix_free(matrix_inverse2);
    matrix_t *matrix_inverse = matrix_solve_cholesky_f(matrix, Id);
    matrix_free(Id);
    return(matrix_inverse);
}

double matrix_det_cholesky_f(const matrix_t *matrix)
{
    if(!sanity_check((void *)matrix, __func__))return 0;
    if(!square_check(matrix, __func__))return 0; 
    double complex *L = matrix_cholesky_f(matrix);
    complex det = 1;
    int n = matrix->rows;
    for (int i=0; i < n; i++)det *= L[i*n+i];
    free(L);
    return det * det;
}