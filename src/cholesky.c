#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "matrix.h"
#include "tools.h"
#include "check.h"

static TYPE complex * matrix_cholesky_f(const matrix_t *matrix);


static TYPE complex * matrix_cholesky_f(const matrix_t *matrix)
{
    if(!sanity_check((void *)matrix, __func__))return NULL;
    if(!square_check(matrix, __func__))return NULL;
    if(!symetry_check(matrix, __func__))return NULL;
    int n = matrix->rows;
    TYPE complex sum;
    size_t size = n*n*sizeof(TYPE complex);
    TYPE complex *L = aligned_alloc(32, size);
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

TYPE complex ** matrix_solve_diag_inf_comp(TYPE complex *A[], TYPE complex *B[], int n, int m)
{
    TYPE complex **X = calloc(n,sizeof(TYPE complex *));
    if(!X){
        perror(__func__);
        return NULL;
    }
    for (int i = 0; i < n; i++){
        X[i] = calloc(m,sizeof(TYPE complex));
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

TYPE complex ** matrix_solve_diag_sup_comp(TYPE complex *A[], TYPE complex *B[], int n, int m)
{
    TYPE complex **X = calloc(n,sizeof(TYPE complex *));
    if(!X){
        perror(__func__);
        return NULL;
    }
    for (int i = 0; i < n; i++){
        X[i] = calloc(m,sizeof(TYPE complex));
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
    TYPE complex **comp_B = calloc(n,sizeof(TYPE complex *));
    if(!comp_B){
        perror(__func__);
        return NULL;
    }
    for (int i = 0; i < n; i++){
        comp_B[i] = calloc(m,sizeof(TYPE complex));
    }
    for (int i=0; i < n; i++){
        #pragma omp simd
        for (int j=0; j < m; j++){
            comp_B[i][j] = (TYPE complex)B->coeff[i][j];
        }
    }
    TYPE complex *L = matrix_cholesky_f(A);
    TYPE complex **stackL = calloc(n,sizeof(TYPE complex *));
    if(!stackL){
        perror(__func__);
        return NULL;
    }
    TYPE complex **stackLT = calloc(n,sizeof(TYPE complex *));
    if(!stackLT){
        perror(__func__);
        return NULL;
    }
    for (int i = 0; i < n; i++){
        stackL[i] = calloc(n,sizeof(TYPE complex));
        stackLT[i] = calloc(n,sizeof(TYPE complex));
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            stackL[i][j] = L[i*n+j];
            stackLT[j][i] = L[i*n+j];
        }
    }
    free(L);
    TYPE complex **Z = matrix_solve_diag_inf_comp(stackL, comp_B, n, m);
    for (int i = 0; i < n; i++){
        free(comp_B[i]);
        free(stackL[i]);
    }
    free(comp_B);
    free(stackL);
    TYPE complex **X = matrix_solve_diag_sup_comp(stackLT, Z, n, m);
    for (int i = 0; i < n; i++){
        free(Z[i]);
        free(stackLT[i]);
    }
    free(Z);
    free(stackLT);
    matrix_t *ret = matrix_create(n,m);
    if(ret){
        for (int i=0; i < n; i++){
            #pragma omp simd
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
    TYPE complex *L = matrix_cholesky_f(matrix);
    complex det = 1;
    int n = matrix->rows;
    for (int i=0; i < n; i++)
        det *= L[i*n+i];
    free(L);
    return det * det;
}
