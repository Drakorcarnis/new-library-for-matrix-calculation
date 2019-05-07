#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "tools.h"
#include "check.h"
// Methods based upon raw determinant calculation. For fun only. Do never use them, cuz you've NO reason to use them. Really.
static matrix_t * matrix_shrink_f(const matrix_t *matrix, unsigned int skipped_row, unsigned int skipped_column)
{
    if(!sanity_check((void *)matrix, __func__))return NULL;
    unsigned int l = 0;
    matrix_t *shrink_matrix = matrix_create(matrix->rows-1, matrix->columns-1);
    for (unsigned int i = 0; i < shrink_matrix->rows; i++) {
        if (l == skipped_row) l++;
        unsigned int k=0;
        for (unsigned int j = 0; j < shrink_matrix->columns; j++) {
            if (k == skipped_column) k++;
            shrink_matrix->coeff[i][j] = matrix->coeff[l][k++];
        }
        l++;
    }
    return shrink_matrix;
}
 
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
