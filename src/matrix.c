#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "matrix.h"
#include "tools.h"
#include "check.h"

// Matrix creation functions

matrix_t * matrix_create(size_t rows, size_t columns)
{
    size_t i = 0;
    matrix_t *matrix = malloc(sizeof(matrix_t));
    if (!matrix)
        goto failed_matrix;
    matrix->rows = rows;
    matrix->columns = columns;
    matrix->coeff = malloc(rows*sizeof(TYPE *));
    if (!matrix->coeff)
        goto failed_coeff;
    size_t size = (columns + 4-columns%4)* sizeof(TYPE);
    for (i=0; i < rows; i++){
        matrix->coeff[i] = aligned_alloc(64, size);
        if (!matrix->coeff[i])
            goto failed_coeff_elt;
        memset(matrix->coeff[i], 0, size);
    }
    return matrix;
failed_coeff_elt:
    for (size_t j = 0; j < i; j++)
        free(matrix->coeff[j]);
    free(matrix->coeff);
failed_coeff:
    free(matrix);
failed_matrix:
    perror(__func__);
    return NULL;
}

matrix_t * matrix_identity(size_t n)
{
    matrix_t * matrix = matrix_create(n, n);
    if(matrix){
        for (size_t i = 0; i < n; i++)
            matrix->coeff[i][i] = 1;
    }
    return matrix;
}

matrix_t * matrix_copy(const matrix_t *matrix)
{
    matrix_t *copy = matrix_create(matrix->rows, matrix->columns);
    if(copy){
        for (size_t i = 0; i < matrix->rows; i++) {
            #pragma omp simd
            for (size_t j = 0; j < matrix->columns; j++) {
                copy->coeff[i][j] = matrix->coeff[i][j];
            }
        }
    }
    return copy;
}

void matrix_free(matrix_t *matrix)
{
    if(!sanity_check(matrix, __func__))return; 
    for (size_t i = 0; i < matrix->rows; i++) {
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
    if(transpose_matrix){
        #pragma omp parallel for
        for (size_t i = 0; i < transpose_matrix->rows; i++)
            for (size_t j = 0; j < transpose_matrix->columns; j++)
                transpose_matrix->coeff[i][j]=matrix->coeff[j][i];
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
    if(add_matrix){
        for (size_t i = 0; i < add_matrix->rows; i++) {
            #pragma omp simd
            for (size_t j = 0; j < add_matrix->columns; j++) {
                add_matrix->coeff[i][j] = matrix1->coeff[i][j] + matrix2->coeff[i][j];
            }
        }
    }
    return add_matrix;
}

matrix_t * matrix_mult_scalar_f(const matrix_t *matrix, TYPE lambda)
{
    if(!sanity_check((void *)matrix, __func__))return NULL;
    matrix_t *mult_matrix = matrix_create(matrix->rows, matrix->columns);
    if(mult_matrix){
        for (size_t i = 0; i < matrix->rows; i++) {
            #pragma omp simd
            for (size_t j = 0; j < matrix->columns; j++) {
                mult_matrix->coeff[i][j]=lambda * matrix->coeff[i][j];
            }
        }
    }
    return mult_matrix;
}
    
matrix_t * matrix_mult_f(const matrix_t *matrix1, const matrix_t *matrix2)
{
    if(!sanity_check((void *)matrix1, __func__))return NULL; 
    if(!sanity_check((void *)matrix2, __func__))return NULL; 
    if((matrix2->rows != matrix1->columns)){
        fprintf(stderr, "%s: not multiplicable matrix (matrix2->rows != matrix1->columns)\n", __func__);
        return NULL;
    }
   	size_t n = matrix1->rows, m = matrix2->columns;
    matrix_t *mult = matrix_create(matrix1->columns, matrix2->columns);
    if(!sanity_check((void *)mult, __func__))return NULL; 
    matrix_t *columns = matrix_transp_f(matrix2);
    if(!sanity_check((void *)columns, __func__))return NULL; 
    size_t default_step = 2*sysconf(_SC_LEVEL1_DCACHE_LINESIZE)/sizeof(TYPE);
    size_t step = default_step > 0 ? default_step:16;
    #pragma omp parallel for
    for (size_t i = 0; i < n; i+=step) {
        int ie = n < i+step ? n : i+step;
        for (size_t j = 0; j < m; j+=step) {
            int je = m < j+step ? m : j+step;
            for (int ii = i; ii < ie; ii++){
                for (int jj = j; jj < je; jj++){
                    TYPE sum = 0;
                    for (size_t k = 0; k < n; k++)
                        sum+=matrix1->coeff[ii][k]*columns->coeff[jj][k];
                    mult->coeff[ii][jj] += sum;
                }
            }
        }
    }
    matrix_free(columns);
    return mult;
}
// matrix_t * matrix_mult_f(const matrix_t *matrix1, const matrix_t *matrix2)
// {
// ... Initialize mul1 and mul2
// int i, i2, j, j2, k, k2;
// double *restrict rres;
// double *restrict rmul1;
// double *restrict rmul2;
    // size_t default_step = 2*sysconf(_SC_LEVEL1_DCACHE_LINESIZE)/sizeof(TYPE);
    // size_t step = default_step > 0 ? default_step:16;
// for (i = 0; i < N; i += SM)
    // for (j = 0; j < N; j += SM)
        // for (k = 0; k < N; k += SM)
            // for (i2 = 0, rres = &res[i][j], rmul1 = &mul1[i][k]; i2 < SM;++i2, rres += N, rmul1 += N){
                // _mm_prefetch (&rmul1[8], _MM_HINT_NTA);
                // for (k2 = 0, rmul2 = &mul2[k][j]; k2 < SM; ++k2, rmul2 += N){
                    // __m128d m1d = _mm_load_sd (&rmul1[k2]);
                    // m1d = _mm_unpacklo_pd (m1d, m1d);
                    // for (j2 = 0; j2 < SM; j2 += 2){
                        // __m128d m2 = _mm_load_pd (&rmul2[j2]);
                        // __m128d r2 = _mm_load_pd (&rres[j2]);
                        // _mm_store_pd (&rres[j2],
                        // _mm_add_pd (_mm_mul_pd (m2, m1d), r2));
                    // }
                // }
            // }
// ... use res matrix
// return 0;
// }


matrix_t * matrix_pow_f(const matrix_t *matrix, int pow)
{
    if(!sanity_check((void *)matrix, __func__))return NULL;  
    if(!square_check(matrix, __func__))return NULL;
    matrix_t *tmp_matrix0 = matrix_copy(matrix);
    matrix_t *tmp_matrix1;
    if(tmp_matrix0){
        for (int i = 0; i < pow-1; i++) {
            tmp_matrix1 = matrix_mult_f(tmp_matrix0, matrix);
            if(tmp_matrix1){
                matrix_free(tmp_matrix0);
                tmp_matrix0 = tmp_matrix1;
            }
        }
    }
    return tmp_matrix0;
}
