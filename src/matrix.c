#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/sysinfo.h>
#include "matrix.h"
#include "tools.h"
#include "check.h"
#include "thread_pool.h"

// Matrix creation functions
thread_pool_t thread_pool;
int libmatrix_init(void)
{
    //Creating thread pool
    int nthreads = 3*get_nprocs()/2;
    // int nthreads = 2;
    printf("Création d'un pool de \x1b[36m%d\x1b[0m threads\n", nthreads);
    if(thread_pool_create(&thread_pool, nthreads, NULL) != THREAD_POOL_OK){
        printf("\x1b[31mproblem0\x1b[0m\n");
        return 0;
    }
    return 1;
}

int libmatrix_end(void)
{
    //Creating thread pool
    printf("Destruction du pool de threads\n");
    if(thread_pool_destroy(&thread_pool) != THREAD_POOL_OK){
        printf("\x1b[31mproblem1\x1b[0m\n");
        return 0;
    }
    return 1;
}

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
        matrix->coeff[i] = aligned_alloc(ALIGN, size);
        // if(posix_memalign((void *)&matrix->coeff[i], ALIGN, size) != 0){
            // perror(NULL);
            // return NULL;
        // }
        // matrix->coeff[i] = malloc(size);
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
    if(!sanity_check((void *)matrix, __func__))return NULL; 
    matrix_t *copy = matrix_create(matrix->rows, matrix->columns);
    if(copy){
        for (size_t i = 0; i < matrix->rows; i++) {
            copy->coeff[i] = __builtin_assume_aligned(copy->coeff[i], ALIGN);
            matrix->coeff[i] = __builtin_assume_aligned(matrix->coeff[i], ALIGN);
            #pragma GCC ivdep
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
typedef struct {
    matrix_t *transpose_matrix;
    const matrix_t *matrix;
}transpose_arg_t;

static void _transpose_task(void* args, int i)
{
    transpose_arg_t *arg = args;
    matrix_t *transpose_matrix = arg->transpose_matrix;
    const matrix_t *matrix = arg->matrix;
    transpose_matrix->coeff[i] = __builtin_assume_aligned(transpose_matrix->coeff[i], ALIGN);
    matrix->coeff[i] = __builtin_assume_aligned(matrix->coeff[i], ALIGN);
    for (size_t j = 0; j < transpose_matrix->columns; j++)
       transpose_matrix->coeff[i][j] = matrix->coeff[j][i];
}

matrix_t * matrix_transp_f(const matrix_t *matrix)
{
    if(!sanity_check((void *)matrix, __func__))return NULL; 
    matrix_t *transpose_matrix = matrix_create(matrix->columns, matrix->rows);
    if(!transpose_matrix)return NULL;
    transpose_arg_t arg = {transpose_matrix, matrix};
    thread_pool_work_t work = {0, NULL, _transpose_task, (void *)&arg};
    for (size_t i = 0; i < transpose_matrix->rows; i++){
        thread_pool_queue_work(&thread_pool, &work, i);
    }
    if(thread_pool_wait(&thread_pool) != THREAD_POOL_OK){
        printf("\x1b[31mproblem\x1b[0m\n");
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
            add_matrix->coeff[i] = __builtin_assume_aligned(add_matrix->coeff[i], ALIGN);
            matrix1->coeff[i] = __builtin_assume_aligned(matrix1->coeff[i], ALIGN);
            matrix2->coeff[i] = __builtin_assume_aligned(matrix2->coeff[i], ALIGN);
            #pragma GCC ivdep
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
            mult_matrix->coeff[i] = __builtin_assume_aligned(mult_matrix->coeff[i], ALIGN);
            matrix->coeff[i] = __builtin_assume_aligned(matrix->coeff[i], ALIGN);
            #pragma GCC ivdep
            for (size_t j = 0; j < matrix->columns; j++) {
                mult_matrix->coeff[i][j]=lambda * matrix->coeff[i][j];
            }
        }
    }
    return mult_matrix;
}

typedef struct {
    size_t n, m, step;
    const matrix_t *matrix1;
    matrix_t *mult, *columns;
}mult_work_t;

static void _mult_task(void* arg, int i)
{
    mult_work_t *work = arg;
    size_t n = work->n;
    size_t m = work->m;
    size_t step = work->step;
    TYPE **matrix_coeff = work->matrix1->coeff;
    TYPE **mult_coeff = work->mult->coeff;
    TYPE **columns_coeff = work->columns->coeff;
    int ie = n < i+step ? n : i+step;
    for (size_t j = 0; j < m; j+=step) {
        int je = m < j+step ? m : j+step;
        for (int ii = i; ii < ie; ii++){
            matrix_coeff[ii] = __builtin_assume_aligned(matrix_coeff[ii], ALIGN);
            for (int jj = j; jj < je; jj++){
                TYPE sum = 0;
                columns_coeff[jj] = __builtin_assume_aligned(columns_coeff[jj], ALIGN);
                for (size_t k = 0; k < n; k++)
                    sum+=matrix_coeff[ii][k]*columns_coeff[jj][k];
                mult_coeff[ii][jj] += sum;
            }
        }
    }
}
matrix_t * matrix_mult_f(const matrix_t *matrix1, const matrix_t *matrix2)
{
    if(!sanity_check((void *)matrix1, __func__))return NULL; 
    if(!sanity_check((void *)matrix2, __func__))return NULL; 
    if((matrix2->rows != matrix1->columns)){
        fprintf(stderr, "%s: not multiplicable matrix (matrix2->rows != matrix1->columns)\n", __func__);
        return NULL;
    }
    matrix_t *mult = matrix_create(matrix1->columns, matrix2->columns);
    if(!sanity_check((void *)mult, __func__))return NULL;
    matrix_t *columns = matrix_transp_f(matrix2);
    if(!sanity_check((void *)columns, __func__))return NULL; 
    size_t default_step = 2*sysconf(_SC_LEVEL1_DCACHE_LINESIZE)/sizeof(TYPE);
    size_t step = default_step > 0 ? default_step:16;
    mult_work_t args = {matrix1->rows, matrix2->columns, step, matrix1, mult, columns};
    thread_pool_work_t work = {0, NULL, _mult_task, (void *)&args};
    for (size_t i = 0; i < matrix1->rows; i+=step)
        thread_pool_queue_work(&thread_pool, &work, i);
    if(thread_pool_wait(&thread_pool) != THREAD_POOL_OK){
        printf("\x1b[31mproblem\x1b[0m\n");
    }
    matrix_free(columns);
    return mult;
}

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
