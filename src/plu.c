#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "matrix.h"
#include "tools.h"
#include "check.h"

typedef struct {
    size_t nb_perm;
    size_t (*perm)[2];
    matrix_t *L;
    matrix_t *U;
} plu_t;

static plu_t * plu_create(size_t rank);
static void plu_free(plu_t *plu);
static plu_t * matrix_plu_f(const matrix_t *matrix);

static plu_t * plu_create(size_t rank){
    plu_t *plu = calloc(1, sizeof(plu_t));
    if(!plu){
        perror(__func__);
        return NULL;
    }
    plu->L = matrix_create(rank, rank);
    plu->U = matrix_identity(rank);
    plu->nb_perm = 0;
    return plu;
}

static void plu_free(plu_t *plu)
{
    if(!sanity_check(plu, __func__))return; 
    matrix_free(plu->L);
    matrix_free(plu->U);
    if(plu->perm)
        free(plu->perm);
    free(plu);   
}

static inline void matrix_row_permute(matrix_t *matrix, int i, int j)
{
    TYPE *tmp = matrix->coeff[i];
    matrix->coeff[i] = matrix->coeff[j];
    matrix->coeff[j] = tmp;
}

static plu_t * matrix_plu_f(const matrix_t *matrix)
{
    if(!sanity_check((void *)matrix, __func__))return NULL;
    if(!square_check(matrix, __func__))return NULL;
    size_t p, n = matrix->rows;
    plu_t *plu = plu_create(n); 
    if(!sanity_check((void *)plu, __func__))return NULL;
    matrix_t *M = matrix_copy(matrix);
    if(!sanity_check((void *)M, __func__))return NULL;
    TYPE **A = M->coeff, **L = plu->L->coeff, **U = plu->U->coeff;
    TYPE sum;
    size_t default_step = 2*sysconf(_SC_LEVEL1_DCACHE_LINESIZE)/sizeof(TYPE);
    size_t step = default_step > 0 ? default_step:16;
    long long time = mstime();
    for (size_t i = 0; i < n; i++){
        for (size_t j = i; j < n; j+=step){
            int je = n < j+step ? n : j+step;
            for (int jj = j; jj < je; jj++){
                sum = 0;
                for (size_t k = 0; k < i; k++)
                    sum += L[jj][k] * U[i][k];
                L[jj][i] = A[jj][i] - sum;
            }
        }
        p = i;
        while((p < n-1) && L[p][i] == 0 && ++p);
        if (p != i){
            matrix_row_permute(plu->L, p, i);
            matrix_row_permute(M, p, i);
            plu->nb_perm++;
            plu->perm = plu->nb_perm > 1 ? realloc(plu->perm,plu->nb_perm*sizeof(*plu->perm)):malloc(plu->nb_perm*sizeof(*plu->perm));
            plu->perm[plu->nb_perm-1][0] = p;
            plu->perm[plu->nb_perm-1][1] = i;
        }
        for (size_t j = i+1; j < n; j+=step){
            int je = n < j+step ? n : j+step;
            for (int jj = j; jj < je; jj++){
                sum = 0;
                for (size_t k = 0; k < i; k++)
                    sum += L[i][k] * U[jj][k];
                U[jj][i] = (A[i][jj] - sum) / L[i][i];
            }
        }
    }
    printf("plu: %s\n", format_time(mstime()-time, "ms"));
    matrix_free(M);
    matrix_t *trueU = matrix_transp_f(plu->U);
    matrix_free(plu->U);
    plu->U = trueU;
    return(plu);  
}

static matrix_t * matrix_solve_low_trig(const matrix_t *A, const matrix_t *B)
{
    // La ruse du B en ligne et chaque element permet de calculer 1 ligne
    if(!sanity_check((void *)A, __func__))return NULL;
    if(!square_check(A, __func__))return NULL; 
    size_t n = A->rows;
    size_t m = B->columns;
    size_t default_step = 2*sysconf(_SC_LEVEL1_DCACHE_LINESIZE)/sizeof(TYPE);
    size_t step = default_step > 0 ? default_step:16;
    matrix_t *X = matrix_transp_f(B);
    for (size_t i = 0; i < m; i+=step){
        int ie = m < i+step ? m : i+step;
        for (size_t j = 0; j < n; j+=step){
            int je = n < j+step ? n : j+step;
            for (int ii = i; ii < ie; ii++){
                for (int jj = j; jj < je; jj++){
                    TYPE sum = X->coeff[ii][jj];
                    for (int k = 0; k < jj; k++)
                        sum -= X->coeff[ii][k] * A->coeff[jj][k];
                    X->coeff[ii][jj] = sum / A->coeff[jj][jj];
                }
            }
        }
    }
    return(matrix_transp_f(X));
}

static matrix_t * matrix_solve_up_trig(const matrix_t *A, const matrix_t *B)
{
    if(!sanity_check((void *)A, __func__))return NULL;
    if(!square_check(A, __func__))return NULL; 
    size_t n = A->rows;
    size_t m = B->columns;
    size_t default_step = 2*sysconf(_SC_LEVEL1_DCACHE_LINESIZE)/sizeof(TYPE);
    size_t step = default_step > 0 ? default_step:16;
    matrix_t *X = matrix_transp_f(B);
    for (size_t i = 0; i < m; i+=step){
        int ie = m < i+step ? m : i+step;
        for (int j = n - 1; j >= 0; j-=step){
            int je = 0 >= j-(int)step ? 0 : j-step;
            for (int ii = i; ii < ie; ii++){
                for (int jj = j-1; jj >= je; jj--){
                    TYPE sum = X->coeff[ii][jj];
                    for (int k = n - 1; k > jj; k--)
                        sum -= X->coeff[ii][k] * A->coeff[jj][k];
                    X->coeff[ii][jj] = sum / A->coeff[jj][jj];
                }
            }
        }
    }
    return(matrix_transp_f(X));
}
TYPE matrix_det_plu_f(const matrix_t *matrix)
{
    if(!sanity_check((void *)matrix, __func__))return 0;
    if(!square_check(matrix, __func__))return 0; 
    plu_t *plu = matrix_plu_f(matrix);
    TYPE det = pow(-1.0, plu->nb_perm);
    for (size_t i=0; i < plu->L->rows; i++)
        det *= plu->L->coeff[i][i];
    plu_free(plu);
    return det;
}
  
matrix_t * matrix_solve_plu_f(const matrix_t *A, const matrix_t *B)
{
    if(!sanity_check((void *)A, __func__))return NULL;
    if(!square_check(A, __func__))return NULL; 
    plu_t *plu = matrix_plu_f(A);
    long long time = mstime();
    matrix_t *permB = matrix_copy(B);
    for (size_t i = 0; i < plu->nb_perm; i++)
        matrix_row_permute(permB, plu->perm[i][0], plu->perm[i][1]);
    printf("perm: %s\n", format_time(mstime()-time, "ms"));
    time = mstime();
    matrix_t *Z = matrix_solve_low_trig(plu->L, permB);
    printf("diag inf: %s\n", format_time(mstime()-time, "ms"));
    matrix_free(permB);
    time = mstime();
    matrix_t *X = matrix_solve_up_trig(plu->U, Z);
    printf("diag sup: %s\n", format_time(mstime()-time, "ms"));
    matrix_free(Z);
    plu_free(plu);
    return(X);
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
