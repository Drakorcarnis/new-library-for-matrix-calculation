#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <immintrin.h>
#include <omp.h>
#include "matrix.h"
#include "matrix_tools.h"
#include "check.h"



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
    perror(__func__);
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
    if(!sanity_check((void *)mult, __func__))return NULL; 
    matrix_t *columns = matrix_transp_f(matrix2);
    if(!sanity_check((void *)columns, __func__))return NULL; 
    double ** coeff = mult->coeff;
    #pragma omp parallel private(i, j) num_threads(8)
	{
        #pragma omp for schedule(static)
        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                coeff[i][j] = avx_mult(n, matrix1->coeff[i], columns->coeff[j]);
            }
        }
    }
    matrix_free(columns);
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


