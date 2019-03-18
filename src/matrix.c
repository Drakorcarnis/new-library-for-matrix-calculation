#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "includes/matrix.h"

static int sanity_check(const void *pointer, const char * function_name)
{
	if(!pointer)
    {
		fprintf(stderr, "%s: NULL pointer\n",function_name);
		return 0;
	}
	return 1;
}

static int square_check(const matrix_t *matrix, const char * function_name)
{
	if((matrix->rows != matrix->columns))
    {
		fprintf(stderr, "%s: not square matrix\n",function_name);
		return 0;
	}
	return 1;
}

static double str2double(const char *str, int len)
{
    char *s = calloc(len+1, sizeof(char));
    memcpy(s, str, len);
    double ret = strtod(s,NULL);
    free(s);
    return ret;
}

matrix_t * matrix_create(int rows, int columns)
{
    int i,j;
    matrix_t *matrix = malloc(sizeof(matrix_t));
    matrix->rows = rows;
    matrix->columns = columns;
    matrix->coeff = malloc(rows*sizeof(double *));
    for (i = 0; i < rows; i++) {
        matrix->coeff[i] = malloc(columns*sizeof(double));
        if(!sanity_check((void *)matrix->coeff[i], __func__))return NULL;
		for (j = 0; j < columns; j++) {
			matrix->coeff[i][j] = 0;
		}
    }
    return matrix;
}

void matrix_free(matrix_t *matrix)
{
    if(!sanity_check(matrix, __func__))return; 
    int i;
    for (i = 0; i < matrix->rows; i++) {
        free(matrix->coeff[i]);
    }
    free(matrix->coeff);
    free(matrix);   
}

plu_t * plu_create(int rank){
    plu_t *plu = malloc(sizeof(plu_t));
    plu->P = matrix_identity(rank);
    plu->L = matrix_create(rank, rank);
    plu->U = matrix_identity(rank);
    plu->nb_perm = 0;
    return plu;
}

void plu_free(plu_t *plu)
{
    if(!sanity_check(plu, __func__))return; 
    matrix_free(plu->P);
    matrix_free(plu->L);
    matrix_free(plu->U);
    free(plu);   
}

matrix_t * matrix_copy(const matrix_t *matrix)
{
    matrix_t *copy = matrix_create(matrix->rows, matrix->columns);
    for (int i = 0; i < matrix->rows; i++) {
		for (int j = 0; j < matrix->columns; j++) {
			copy->coeff[i][j] = matrix->coeff[i][j];
		}
    }
    return copy;
}


matrix_t * str2matrix(int argc, char **argv, char separator)
{
    int i, j, columns = 0, len, trigger, flag = 0, count;
    matrix_t *matrix = NULL;
    for (i = 0; i < argc; i++) {
        count = 0;
        for (j = 0; argv[i][j]; j++){
            if(argv[i][j] != separator && argv[i][j] != '\t' && argv[i][j] != '\n'){
                if(!flag)count++;
                flag = 1;
            } else {
                flag = 0;
            }
        }
        if (count>columns)columns=count;
    }
    matrix = matrix_create(argc, columns);
    if(!sanity_check((void *)matrix, __func__))return NULL; 
    for (i = 0; i < matrix->rows; i++) {
        len=0;
        flag=0;
        trigger=0;
        for (j = 0; j <= (int)strlen(argv[i]); j++) {
            if(argv[i][j] != separator && argv[i][j] != '\t' && argv[i][j] != '\n'){
                if(!flag)trigger=j;
                flag = 1;
            } else {
                if(flag)matrix->coeff[i][len++]=str2double(argv[i]+trigger,j-trigger);
                flag = 0;
            }
        }
        if(trigger < j && flag && len < columns)matrix->coeff[i][len]=str2double(argv[i]+trigger,j-trigger-1);
    }
    return matrix;
}

matrix_t * file2matrix(char * filename)
{
    FILE *fp = fopen(filename, "r");
	int argc = 0, i;
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	if(!sanity_check(fp, __func__))return NULL;
	char **argv = malloc(sizeof(char*));
	while ((read = getline(&line, &len, fp)) != -1) {
		if (read > 1 || ((read == 1) && (strcmp(line, "\n") != 0))){
            argv = realloc(argv, (argc + 1) * sizeof(char*));
			argv[argc] = calloc(read+1, sizeof(char));
			strncpy(argv[argc++], line, read);
		}	
    }
    if(line)free(line);
    fclose(fp);
	matrix_t *matrix = str2matrix(argc, argv, ' ');
	for (i=0; i<argc; i++)
	{
		free(argv[i]);
	}
	free(argv);
	return matrix;
}

void matrix_display(const matrix_t *matrix)
{
    if(!sanity_check((void *)matrix, __func__))return; 
    int i, j;
    for (i = 0; i < matrix->rows; i++) {
        printf("[");
        for (j = 0; j < matrix->columns; j++){
            if(matrix->coeff[i][j] < 1e-10)
                printf("0 ");
            else
                printf("%g ", matrix->coeff[i][j]);
    }
        printf("]\n");
    }
}

matrix_t * matrix_transp_f(const matrix_t *matrix)
{
    if(!sanity_check((void *)matrix, __func__))return NULL; 
    int i, j;
    matrix_t *transpose_matrix = matrix_create(matrix->columns, matrix->rows);
    for (i = 0; i < transpose_matrix->rows; i++) {
        for (j = 0; j < transpose_matrix->columns; j++)transpose_matrix->coeff[i][j]=matrix->coeff[j][i];
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
    int i, j;
    matrix_t *add_matrix = matrix_create(matrix1->rows, matrix1->columns);
    for (i = 0; i < add_matrix->rows; i++) {
        for (j = 0; j < add_matrix->columns; j++) {
                add_matrix->coeff[i][j] = matrix1->coeff[i][j] + matrix2->coeff[i][j];
        }
    }
    return add_matrix;
}

matrix_t * matrix_mult_scalar_f(const matrix_t *matrix, double lambda)
{
    if(!sanity_check((void *)matrix, __func__))return NULL;
    int i, j;
    matrix_t *mult_matrix = matrix_create(matrix->rows, matrix->columns);
    for (i = 0; i < matrix->rows; i++) {
        for (j = 0; j < matrix->columns; j++) {
                mult_matrix->coeff[i][j]=lambda * matrix->coeff[i][j];
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
	int i, j, k;
    matrix_t *mult = matrix_create(matrix1->columns, matrix2->columns);
    for (i = 0; i < mult->rows; i++) {
        for (j = 0; j < mult->columns; j++) {
            for (k = 0; k < matrix1->rows; k++) {
                mult->coeff[i][j] = mult->coeff[i][j] + (matrix1->coeff[i][k] * matrix2->coeff[k][j]);
            }
        }
    }
    return mult;
}

matrix_t * matrix_pow_f(const matrix_t *matrix, int pow)
{
    if(!sanity_check((void *)matrix, __func__))return NULL;  
	if(!square_check(matrix, __func__))return NULL; 
	int i, j, k;
    matrix_t *pow_matrix = matrix_create(matrix->rows, matrix->columns);
    matrix_t *tmp_matrix;
	for (j = 0; j < pow_matrix->rows; j++) {
		for (k = 0; k < pow_matrix->columns; k++) {
			pow_matrix->coeff[j][k] = matrix->coeff[j][k];
		}
	}
    for (i = 0; i < pow-1; i++) {
		tmp_matrix = matrix_mult_f(pow_matrix, pow_matrix);
		for (j = 0; j < pow_matrix->rows; j++) {
			for (k = 0; k < pow_matrix->columns; k++) {
				pow_matrix->coeff[j][k] = tmp_matrix->coeff[j][k];
			}
		}
		matrix_free(tmp_matrix);
    }
    return pow_matrix;
}


matrix_t * matrix_shrink_f(const matrix_t *matrix, int skipped_row, int skipped_column)
{
    if(!sanity_check((void *)matrix, __func__))return NULL;
    int i, j, k, l=0;
    matrix_t *shrink_matrix = matrix_create(matrix->rows-1, matrix->columns-1);
    for (i = 0; i < shrink_matrix->rows; i++) {
		if (l == skipped_row) l++;
        k=0;
        for (j = 0; j < shrink_matrix->columns; j++) {
            if (k == skipped_column) k++;
            shrink_matrix->coeff[i][j]=matrix->coeff[l][k++];
        }
		l++;
    }
    return shrink_matrix;
}

double matrix_det_raw_f(const matrix_t *matrix)
{
    int i;
    double det = 0;
    double sign = 0;
    matrix_t *shrinked_matrix = NULL;
    if(!sanity_check((void *)matrix, __func__))return 0;
    if(matrix->columns != matrix->rows) return 0;
    if(matrix->columns == 1)return matrix->coeff[0][0];
    for (i = 0; i < matrix->columns; i++) 
    {
        shrinked_matrix=matrix_shrink_f(matrix, 0, i);
		sign = (int)pow(-1.0,(double)i);
        det += sign * matrix->coeff[0][i] * matrix_det_raw_f(shrinked_matrix);
        matrix_free(shrinked_matrix);
    }
    return det;
}

double matrix_det_plu_f(const matrix_t *matrix)
{
	if(!sanity_check((void *)matrix, __func__))return 0;
    if(!square_check(matrix, __func__))return 0; 
    plu_t *plu = matrix_plu_f(matrix);
    double det = pow(-1.0, plu->nb_perm);
    for (int i=0; i < plu->L->rows; i++)det *= plu->L->coeff[i][i];
    for (int i=0; i < plu->U->rows; i++)det *= plu->U->coeff[i][i];
    plu_free(plu);
    return det;
}

matrix_t * matrix_identity(int n)
{
    matrix_t * matrix = matrix_create(n, n);
    for (int i = 0; i < n; i++) matrix->coeff[i][i] = 1;
    return matrix;
}

matrix_t * matrix_permutation(int line1, int line2, int n)
{
    double tmp;
    matrix_t * matrix = matrix_identity(n);
    for (int i = 0; i < n; i++)
    {
        tmp = matrix->coeff[line1][i];
        matrix->coeff[line1][i] = matrix->coeff[line2][i];
        matrix->coeff[line2][i] = tmp;
    }
    matrix_t * transpose = matrix_transp_f(matrix);
    matrix_free(matrix);
    return(transpose);
}

plu_t * matrix_plu_f(const matrix_t *matrix)
{
	if(!sanity_check((void *)matrix, __func__))return NULL;
    if(!square_check(matrix, __func__))return NULL;
    int p, n = matrix->rows;
    plu_t *plu = plu_create(n); 
    matrix_t *A = matrix_copy(matrix), *L = plu->L, *U = plu->U;
    matrix_t *perm, *perm_mult;
    double sum, tmp;
    for (int i = 0; i < n; i++)
    {
        for (int j = i; j < n; j++)
        {
            sum = 0;
            for (int k = 0; k < i; k++)sum += L->coeff[j][k] * U->coeff[k][i];
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
            for (int j = 0; j < n; j++)
            {
                tmp = A->coeff[i][j];
                A->coeff[i][j] = A->coeff[p][j];
                A->coeff[p][j] = tmp;
            }
            for (int j = 0; j < n; j++)
            {
                tmp = L->coeff[i][j];
                L->coeff[i][j] = L->coeff[p][j];
                L->coeff[p][j] = tmp;
            }
        } // Fin des permutations
        for (int j = i+1; j < n; j++)
        {
            sum = 0;
            for (int k = 0; k < i; k++)
            {
                sum += L->coeff[i][k] * U->coeff[k][j];
            }  
            U->coeff[i][j] = (A->coeff[i][j] - sum) / L->coeff[i][i];
        }
    }
    matrix_free(A);
    return(plu);  
}

matrix_t * matrix_com_f(const matrix_t *matrix)
{
	if(!sanity_check((void *)matrix, __func__))return NULL;
    if(!square_check(matrix, __func__))return NULL; 
	int i, j;
	double sign = 0;
    matrix_t *shrinked_matrix = NULL;
    matrix_t *co_matrix = matrix_create(matrix->rows, matrix->columns);
    for (i = 0; i < co_matrix->rows; i++)
    {
        for (j = 0; j < co_matrix->columns; j++)
        {
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

matrix_t * matrix_solve_diag_inf(const matrix_t *A, const matrix_t *B)
{
	if(!sanity_check((void *)A, __func__))return NULL;
    if(!square_check(A, __func__))return NULL; 
    int n = A->rows;
    int m = B->columns;
	matrix_t *X = matrix_create(n, n);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            X->coeff[j][i] = B->coeff[j][i];
            for (int k = 0; k < j; k++)X->coeff[j][i] =  X->coeff[j][i] - X->coeff[k][i] * A->coeff[j][k];
            X->coeff[j][i] = X->coeff[j][i] / A->coeff[j][j];
        }
    }
    return(X);
}

matrix_t * matrix_solve_diag_sup(const matrix_t *A, const matrix_t *B)
{
	if(!sanity_check((void *)A, __func__))return NULL;
    if(!square_check(A, __func__))return NULL; 
    int n = A->rows;
    int m = B->columns;
	matrix_t *X = matrix_create(n, n);
    for (int i = 0; i < m; i++)
    {
        for (int j = n - 1; j >= 0; j--)
        {
            X->coeff[j][i] = B->coeff[j][i];
            for (int k = n - 1; k > j; k--)X->coeff[j][i] =  X->coeff[j][i] - X->coeff[k][i] * A->coeff[j][k];
            X->coeff[j][i] = X->coeff[j][i] / A->coeff[j][j];
        }
    }
    return(X);
}
    
matrix_t * matrix_inverse_plu_f(const matrix_t *matrix)
{
	if(!sanity_check((void *)matrix, __func__))return NULL;
    if(!square_check(matrix, __func__))return NULL; 
	matrix_t *I = matrix_identity(matrix->rows);
    matrix_t *matrix_inverse = matrix_solve_plu_f(matrix, I);
	matrix_free(I);
	return(matrix_inverse);
}  
  
matrix_t * matrix_solve_plu_f(const matrix_t *A, const matrix_t *B)
{
	if(!sanity_check((void *)A, __func__))return NULL;
    if(!square_check(A, __func__))return NULL; 
    plu_t *plu = matrix_plu_f(A);
    matrix_t *Z = matrix_solve_diag_inf(plu->L, B);
    matrix_t *X = matrix_solve_diag_sup(plu->U, Z);
    matrix_t *ret = matrix_mult_f(X, plu->P);
    plu_free(plu);
	matrix_free(X);
	matrix_free(Z);
	return(ret);
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
      