#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "matrix_tools.h"
#include "check.h"

typedef struct {
    int nb_perm;
    matrix_t *P;
    matrix_t *L;
    matrix_t *U;
} plu_t;

static plu_t * plu_create(unsigned int rank);
static void plu_free(plu_t *plu);
static plu_t * matrix_plu_f(const matrix_t *matrix);

static plu_t * plu_create(unsigned int rank){
    plu_t *plu = malloc(sizeof(plu_t));
    if(!plu){
        perror(__func__);
        return NULL;
    }
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
