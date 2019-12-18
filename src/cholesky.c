#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "matrix.h"
#include "tools.h"
#include "check.h"

static COMPLEX_TYPE * matrix_cholesky_f(const matrix_t *matrix);


static COMPLEX_TYPE * matrix_cholesky_f(const matrix_t *matrix)
{
    if(!sanity_check((void *)matrix, __func__))return NULL;
    if(!square_check(matrix, __func__))return NULL;
    if(!symetry_check(matrix, __func__))return NULL;
    int n = matrix->rows;
    COMPLEX_TYPE sum;
    size_t size = n*n*sizeof(COMPLEX_TYPE);
    COMPLEX_TYPE *L = aligned_alloc(ALIGN, size);
    if(!L){
        perror(__func__);
        return NULL;
    }
    memset(L, 0, size);
    for (int i = 0; i < n; i++){
        sum = 0;
        int index = n*i;
        for (int k = 0; k < i; k++)sum += L[index+k] * L[index+k];
        L[index+i] = csqrt(matrix->coeff[i][i] - sum);
        for (int j = i+1; j < n; j++){
            sum = 0;
            int index2 = n*j;
            for (int k = 0; k < i; k++)sum += L[index+k] * L[index2+k];
            L[index2+i] = (matrix->coeff[i][j] - sum)/L[index+i];
        }
    }
    return(L);  
} 

COMPLEX_TYPE ** matrix_solve_diag_inf_comp(COMPLEX_TYPE *A[], COMPLEX_TYPE *B[], int n, int m)
{
    COMPLEX_TYPE **X = calloc(n,sizeof(COMPLEX_TYPE *));
    if(!X){
        perror(__func__);
        return NULL;
    }
    for (int i = 0; i < n; i++){
        X[i] = aligned_alloc(ALIGN, m*sizeof(COMPLEX_TYPE));
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

COMPLEX_TYPE ** matrix_solve_diag_sup_comp(COMPLEX_TYPE *A[], COMPLEX_TYPE *B[], int n, int m)
{
    COMPLEX_TYPE **X = calloc(n,sizeof(COMPLEX_TYPE *));
    if(!X){
        perror(__func__);
        return NULL;
    }
    for (int i = 0; i < n; i++){
        X[i] = aligned_alloc(ALIGN, m*sizeof(COMPLEX_TYPE));
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
    COMPLEX_TYPE **comp_B = calloc(n,sizeof(COMPLEX_TYPE *));
    if(!comp_B){
        perror(__func__);
        return NULL;
    }
    for (int i = 0; i < n; i++){
        comp_B[i] = aligned_alloc(ALIGN, m*sizeof(COMPLEX_TYPE));
    }
    for (int i=0; i < n; i++){
        comp_B[i] = __builtin_assume_aligned(comp_B[i], ALIGN);
        B->coeff[i] = __builtin_assume_aligned(B->coeff[i], ALIGN);
        #pragma GCC ivdep
        for (int j=0; j < m; j++){
            comp_B[i][j] = (COMPLEX_TYPE)B->coeff[i][j];
        }
    }
    COMPLEX_TYPE *L = matrix_cholesky_f(A);
    COMPLEX_TYPE **stackL = calloc(n,sizeof(COMPLEX_TYPE *));
    if(!stackL){
        perror(__func__);
        return NULL;
    }
    COMPLEX_TYPE **stackLT = calloc(n,sizeof(COMPLEX_TYPE *));
    if(!stackLT){
        perror(__func__);
        return NULL;
    }
    for (int i = 0; i < n; i++){
        stackL[i] = aligned_alloc(ALIGN, n*sizeof(COMPLEX_TYPE));
        stackLT[i] = aligned_alloc(ALIGN, n*sizeof(COMPLEX_TYPE));
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            stackL[i][j] = L[i*n+j];
            stackLT[j][i] = L[i*n+j];
        }
    }
    free(L);
    COMPLEX_TYPE **Z = matrix_solve_diag_inf_comp(stackL, comp_B, n, m);
    for (int i = 0; i < n; i++){
        free(comp_B[i]);
        free(stackL[i]);
    }
    free(comp_B);
    free(stackL);
    COMPLEX_TYPE **X = matrix_solve_diag_sup_comp(stackLT, Z, n, m);
    for (int i = 0; i < n; i++){
        free(Z[i]);
        free(stackLT[i]);
    }
    free(Z);
    free(stackLT);
    matrix_t *ret = matrix_create(n,m);
    if(ret){
        for (int i=0; i < n; i++){
            ret->coeff[i] = __builtin_assume_aligned(ret->coeff[i], ALIGN);
            X[i] = __builtin_assume_aligned(X[i], ALIGN);
            #pragma GCC ivdep
            for (int j=0; j < m; j++){
                ret->coeff[i][j] = creal(X[i][j]);
            }
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
    matrix_t *matrix_inverse = matrix_solve_cholesky_f(matrix, Id);
    matrix_free(Id);
    return(matrix_inverse);
}

TYPE matrix_det_cholesky_f(const matrix_t *matrix)
{
    if(!sanity_check((void *)matrix, __func__))return 0;
    if(!square_check(matrix, __func__))return 0; 
    COMPLEX_TYPE *L = matrix_cholesky_f(matrix);
    complex det = 1;
    int n = matrix->rows;
    for (int i=0; i < n; i++)
        det *= L[i*n+i];
    free(L);
    return det * det;
}
