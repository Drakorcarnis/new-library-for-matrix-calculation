#include <stdio.h>
#include "matrix.h"
#include "check.h"

int sanity_check(const void *pointer, const char *function_name)
{
    if(!pointer){
        fprintf(stderr, "%s: NULL pointer\n",function_name);
        return 0;
    }
    return 1;
}

int square_check(const matrix_t *matrix, const char *function_name)
{
    if((matrix->rows != matrix->columns)){
        fprintf(stderr, "%s: not square matrix\n",function_name);
        return 0;
    }
    return 1;
}

int symetry_check(const matrix_t *matrix, const char *function_name)
{
    for (size_t i = 0; i < matrix->rows; i++) {
        for (size_t j = 0; j < i; j++) {
            if(matrix->coeff[i][j] != matrix->coeff[j][i]){
                fprintf(stderr, "%s: not symetric matrix\n",function_name);
                return 0;
            }
        }
    }
    return 1;
}