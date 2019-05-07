#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "matrix.h"
#include "tools.h"

#define DATA_PATH "/home/Drakorcarnis/new_library/src/data"
// Those test must be run with TYPE defined to double in matrix.h
#define PRECISION 11
int precision = PRECISION;

#define IN_1MATRIX_OUT_DOUBLE       matrix_det_plu_f,   matrix_det_cholesky_f,    matrix_det_raw_f
#define IN_1MATRIX_OUT_DOUBLE_IN    "matrix",           "matrix_sym"    ,         "small_matrix"     
#define IN_1MATRIX_OUT_DOUBLE_OUT   41528648144648139187061754753110960955112756453546573882130475948738.000000, 11018591242987731278480209809487334669634103274730489724827165526048.0, -18570.0
#define IN_1MATRIX_OUT_DOUBLE_NAME  "matrix_det_plu_f", "matrix_det_cholesky_f",  "matrix_det_raw_f"

#define IN_1MATRIX_OUT_MATRIX       matrix_transp_f,    matrix_inverse_plu_f,   matrix_inverse_cholesky_f,   matrix_inverse_raw_f,  
#define IN_1MATRIX_OUT_MATRIX_IN    "matrix",           "matrix",               "matrix_sym"  ,              "small_matrix",             
#define IN_1MATRIX_OUT_MATRIX_OUT   "matrix_transp",    "matrix_inv",           "matrix_sym_inv"  ,          "matrix_small_inv",         
#define IN_1MATRIX_OUT_MATRIX_NAME  "matrix_transp_f",  "matrix_inverse_plu_f", "matrix_inverse_cholesky_f", "matrix_inverse_raw_f"

#define IN_2MATRIX_OUT_MATRIX       matrix_add_f,   matrix_mult_f,      matrix_mult_f,      matrix_solve_plu_f,         matrix_solve_cholesky_f,        matrix_solve_plu_f,             matrix_solve_raw_f,  
#define IN_2MATRIX_OUT_MATRIX_IN0   "matrix",       "matrix",           "matrix",           "matrix",                   "matrix_sym",                   "matrix_sym",                   "small_matrix",            
#define IN_2MATRIX_OUT_MATRIX_IN1   "matrix",       "matrix",           "B",                "B",                        "B",                            "B",                            "B_small",    
#define IN_2MATRIX_OUT_MATRIX_OUT   "matrix_add",   "matrix_mult",      "vector_mult",      "X",                        "X_sym",                        "X_sym",                        "X_small",
#define IN_2MATRIX_OUT_MATRIX_NAME  "matrix_add_f", "matrix_mult_f",    "vector_mult_f",    "matrix_solve_plu_f",       "matrix_solve_cholesky_f",      "matrix_sym_solve_plu_f",       "matrix_solve_raw_f",

#define IN_1MATRIX_1DOUBLE_OUT_MATRIX       matrix_mult_scalar_f
#define IN_1MATRIX_1DOUBLE_OUT_MATRIX_IN0   "matrix"
#define IN_1MATRIX_1DOUBLE_OUT_MATRIX_IN1   0.5
#define IN_1MATRIX_1DOUBLE_OUT_MATRIX_OUT   "matrix_mult_scalar"
#define IN_1MATRIX_1DOUBLE_OUT_MATRIX_NAME  "matrix_mult_scalar_f"

#define IN_1MATRIX_1INT_OUT_MATRIX      matrix_pow_f
#define IN_1MATRIX_1INT_OUT_MATRIX_IN0  "matrix"
#define IN_1MATRIX_1INT_OUT_MATRIX_IN1  3
#define IN_1MATRIX_1INT_OUT_MATRIX_OUT  "matrix_pow"
#define IN_1MATRIX_1INT_OUT_MATRIX_NAME "matrix_pow_f"

#define IN_1MATRIX_2INT_OUT_MATRIX      
#define IN_1MATRIX_2INT_OUT_MATRIX_IN0  
#define IN_1MATRIX_2INT_OUT_MATRIX_IN1  
#define IN_1MATRIX_2INT_OUT_MATRIX_IN2  
#define IN_1MATRIX_2INT_OUT_MATRIX_OUT  
#define IN_1MATRIX_2INT_OUT_MATRIX_NAME 

typedef struct {
    char *test_name;
    int result;
    long long mstime;
}result_t;

static void process_result(result_t res)
{
    char *formatted_time = format_time(res.mstime, "ms");
    if (res.result)
        printf("\033[0;32m[OK]");
    else
        printf("\033[0;31m[NOK]");
    printf("%s (%s)\033[0m\n", res.test_name, formatted_time);
    free(formatted_time);
}


static int test_function_in_matrix_out_TYPE_f(TYPE(*function)(const matrix_t*), matrix_t* matrix, TYPE expected, char *test_name)
{
    long long time = mstime();
    TYPE ret = (function)(matrix);
    long long time2 = mstime();
    result_t res = {test_name, fabs(-1+ret/expected) < 2e-3, time2 - time};
    process_result(res);
    if (!res.result){
        printf("input:\n");
        matrix_display_exact(matrix, precision);
        printf("expected: %f\ngot: %f\n",expected, ret);
    }
    return(res.result);
}

static int test_function_in_matrix_out_matrix_f(matrix_t*(*function)(const matrix_t*), matrix_t* matrix, matrix_t* expected, char *test_name)
{
    long long time = mstime();
    matrix_t *ret = (function)(matrix);
    long long time2 = mstime();
    result_t res = {test_name, test_matrix_equality(expected, ret, precision), time2 - time};
    process_result(res);
    if (!res.result){
        printf("input:\n");
        matrix_display_exact(matrix, precision);
        printf("expected:\n");
        matrix_display_exact(expected, precision);
        printf("got:\n");
        matrix_display_exact(ret, precision);
    }
    matrix_free(ret);
    return(res.result);
}

static int test_function_in_2matrix_out_matrix_f(matrix_t*(*function)(const matrix_t*, const matrix_t*), matrix_t* matrix1, matrix_t* matrix2, matrix_t* expected, char *test_name)
{
    long long time = mstime();
    matrix_t *ret = (function)(matrix1, matrix2);
    long long time2 = mstime();
    result_t res = {test_name, test_matrix_equality(expected, ret, precision), time2 - time};
    process_result(res);
    if (!res.result){
        printf("input1:\n");
        matrix_display_exact(matrix1, precision);
        printf("input2:\n");
        matrix_display_exact(matrix2, precision);
        printf("expected:\n");
        matrix_display_exact(expected, precision);
        printf("got:\n");
        matrix_display_exact(ret, precision);
    }
    matrix_free(ret);
    return(res.result);
}

static int test_function_in_matrix_TYPE_out_matrix_f(matrix_t*(*function)(const matrix_t*, TYPE), matrix_t* matrix, TYPE val, matrix_t* expected, char *test_name)
{
    long long time = mstime();
    matrix_t *ret = (function)(matrix, val);
    long long time2 = mstime();
    result_t res = {test_name, test_matrix_equality(expected, ret, precision), time2 - time};
    process_result(res);
    if (!res.result){
        printf("input:\n");
        matrix_display_exact(matrix, precision);
        printf("expected:\n");
        matrix_display_exact(expected, precision);
        printf("got:\n");
        matrix_display_exact(ret, precision);
    }
    matrix_free(ret);
    return(res.result);
}

static int test_function_in_matrix_int_out_matrix_f(matrix_t*(*function)(const matrix_t*, int), matrix_t* matrix, int val, matrix_t* expected, char *test_name)
{
    long long time = mstime();
    matrix_t *ret = (function)(matrix, val);
    long long time2 = mstime();
    result_t res = {test_name, test_matrix_equality(expected, ret, precision), time2 - time};
    process_result(res);
    if (!res.result){
        printf("input:\n");
        matrix_display_exact(matrix, precision);
        printf("expected:\n");
        matrix_display_exact(expected, precision);
        printf("got:\n");
        matrix_display_exact(ret, precision);
    }
    matrix_free(ret);
    return(res.result);
}

static int test_function_in_matrix_2int_out_matrix_f(matrix_t*(*function)(const matrix_t*, size_t, size_t), matrix_t* matrix, int val0, int val1, matrix_t* expected, char *test_name)
{
    long long time = mstime();
    matrix_t *ret = (function)(matrix, val0, val1);
    long long time2 = mstime();
    result_t res = {test_name, test_matrix_equality(expected, ret, precision), time2 - time};
    process_result(res);
    if (!res.result){
        printf("input:\n");
        matrix_display_exact(matrix, precision);
        printf("expected:\n");
        matrix_display_exact(expected, precision);
        printf("got:\n");
        matrix_display_exact(ret, precision);
    }
    matrix_free(ret);
    return(res.result);
}

static int test_error_cases(void)
{
    matrix_t *matrix1 = matrix_create(3,3);
    matrix_t *matrix2 = matrix_create(2,5);
    int res = 1;
    long long time = mstime();
    matrix_t *ret = matrix_add_f(matrix1, matrix2);
    if(ret != NULL){
        matrix_free(ret);
        process_result((result_t){"test_not_addable_matrix", 0, 0});
        res = 0;
    }
    ret = matrix_mult_f(matrix1, matrix2);
    if(ret != NULL){
        matrix_free(ret);
        process_result((result_t){"test_not_multiplicable_matrix", 0, 0});
        res = 0;
    }
    ret = matrix_inverse_raw_f(matrix1);
    if(ret != NULL){
        matrix_free(ret);
        process_result((result_t){"test_not_inversible_matrix", 0, 0});
        res = 0;
    }
    ret = matrix_pow_f(matrix2, 15);
    if(ret != NULL){
        matrix_free(ret);
        process_result((result_t){"test_not_square_matrix", 0, 0});
        res = 0;
    }
    matrix1->coeff[0][1] = 1.0;
    ret = matrix_inverse_cholesky_f(matrix1);
    if(ret != NULL){
        matrix_free(ret);
        process_result((result_t){"test_not_symetric_matrix", 0, 0});
        res = 0;
    }
    ret = matrix_inverse_cholesky_f((matrix_t *)NULL);
    if(ret != NULL){
        matrix_free(ret);
        process_result((result_t){"test_not_sane_matrix", 0, 0});
        res = 0;
    }
    ret = matrix_create(4294967295,4294967295);
    if(ret != NULL){
        free(ret);
        process_result((result_t){"test_failed_alloc", 0, 0});
        res = 0;
    }
    long long time2 = mstime();
    if(res)process_result((result_t){"test_error_cases", 1, time2 - time});
    matrix_free(matrix1);
    matrix_free(matrix2);
    return(res);
}

static matrix_t** chartab2matrixtab(char ** filetab, int size, char *data_path)
{
    matrix_t** matrixtab = malloc(sizeof(matrix_t*) * size);
    for (int i = 0; i < size; i++)
    {
        size_t bufsz = snprintf(NULL, 0, "%s/%s.txt",data_path,filetab[i]);
        char* filename = malloc((bufsz+1)*sizeof(*filename));
        snprintf(filename, bufsz + 1, "%s/%s.txt",data_path,filetab[i]);
        matrixtab[i] = file2matrix(filename);
        free(filename);
    }
    return matrixtab;
}

static void free_matrixtab(matrix_t** matrixtab, int size)
{
    for (int i = 0; i < size; i++)matrix_free(matrixtab[i]);
    free(matrixtab);
    matrixtab = NULL;
}

int main(int argc, char **argv) {
    int ret = 0;
    char *data_path = DATA_PATH;
    if(argc > 1)
        precision = atoi(argv[1]);
    int size;
    
    TYPE(*test_function_in_1matrix_out_TYPE_f[])(const matrix_t*) = {IN_1MATRIX_OUT_DOUBLE};
    char* test_function_in_1matrix_out_TYPE_name[] = {IN_1MATRIX_OUT_DOUBLE_NAME};
    char* test_function_in_1matrix_out_TYPE_in[] = {IN_1MATRIX_OUT_DOUBLE_IN};
    size = sizeof(test_function_in_1matrix_out_TYPE_in)/sizeof(char*);
    matrix_t** test_function_in_1matrix_out_TYPE_in_input = chartab2matrixtab(test_function_in_1matrix_out_TYPE_in, size, data_path);
    TYPE test_function_in_1matrix_out_TYPE_out[] = {IN_1MATRIX_OUT_DOUBLE_OUT}; 
    for (int i = 0; i < size; i++)
        if(!test_function_in_matrix_out_TYPE_f(test_function_in_1matrix_out_TYPE_f[i],test_function_in_1matrix_out_TYPE_in_input[i], test_function_in_1matrix_out_TYPE_out[i], test_function_in_1matrix_out_TYPE_name[i]))ret=1;  
    free_matrixtab(test_function_in_1matrix_out_TYPE_in_input, size);
    
    matrix_t*(*test_function_in_1matrix_out_matrix_f[])(const matrix_t*) = {IN_1MATRIX_OUT_MATRIX};
    char* test_function_in_1matrix_out_matrix_name[] = {IN_1MATRIX_OUT_MATRIX_NAME};
    char* test_function_in_1matrix_out_matrix_in[] = {IN_1MATRIX_OUT_MATRIX_IN};
    size = sizeof(test_function_in_1matrix_out_matrix_in)/sizeof(char*);
    matrix_t** test_function_in_1matrix_out_matrix_in_input = chartab2matrixtab(test_function_in_1matrix_out_matrix_in, size, data_path);
    char* test_function_in_1matrix_out_matrix_out[] = {IN_1MATRIX_OUT_MATRIX_OUT};
    matrix_t** test_function_in_1matrix_out_matrix_out_output = NULL;
    size = sizeof(test_function_in_1matrix_out_matrix_out)/sizeof(char*);
    test_function_in_1matrix_out_matrix_out_output = chartab2matrixtab(test_function_in_1matrix_out_matrix_out, size, data_path);
    for (int i = 0; i < size; i++)
        if(!test_function_in_matrix_out_matrix_f(test_function_in_1matrix_out_matrix_f[i],test_function_in_1matrix_out_matrix_in_input[i], test_function_in_1matrix_out_matrix_out_output[i], test_function_in_1matrix_out_matrix_name[i]))ret=1; 
    free_matrixtab(test_function_in_1matrix_out_matrix_in_input, size);  
    free_matrixtab(test_function_in_1matrix_out_matrix_out_output, size);  
    
    matrix_t*(*test_function_in_2matrix_out_matrix[])(const matrix_t*, const matrix_t*) = {IN_2MATRIX_OUT_MATRIX};
    char* test_function_in_2matrix_out_matrix_name[] = {IN_2MATRIX_OUT_MATRIX_NAME};
    char* test_function_in_2matrix_out_matrix_in0[] = {IN_2MATRIX_OUT_MATRIX_IN0};
    size = sizeof(test_function_in_2matrix_out_matrix_in0)/sizeof(char*);
    matrix_t** test_function_in_2matrix_out_matrix_in0_input = chartab2matrixtab(test_function_in_2matrix_out_matrix_in0, size, data_path);
    char* test_function_in_2matrix_out_matrix_in1[] = {IN_2MATRIX_OUT_MATRIX_IN1};
    size = sizeof(test_function_in_2matrix_out_matrix_in1)/sizeof(char*);
    matrix_t** test_function_in_2matrix_out_matrix_in1_input = chartab2matrixtab(test_function_in_2matrix_out_matrix_in1, size, data_path);
    char* test_function_in_2matrix_out_matrix_out[] = {IN_2MATRIX_OUT_MATRIX_OUT};
    matrix_t** test_function_in_2matrix_out_matrix_out_output = NULL;
    size = sizeof(test_function_in_2matrix_out_matrix_out)/sizeof(char*);
    test_function_in_2matrix_out_matrix_out_output = chartab2matrixtab(test_function_in_2matrix_out_matrix_out, size, data_path);
    for (int i = 0; i < size; i++)
        if(!test_function_in_2matrix_out_matrix_f(test_function_in_2matrix_out_matrix[i],test_function_in_2matrix_out_matrix_in0_input[i], test_function_in_2matrix_out_matrix_in1_input[i], test_function_in_2matrix_out_matrix_out_output[i], test_function_in_2matrix_out_matrix_name[i]))ret=1;
    free_matrixtab(test_function_in_2matrix_out_matrix_in0_input, size);
    free_matrixtab(test_function_in_2matrix_out_matrix_in1_input, size);
    free_matrixtab(test_function_in_2matrix_out_matrix_out_output, size);
    
    matrix_t*(*test_function_in_1matrix_TYPE_out_matrix_f[])(const matrix_t*, TYPE) = {IN_1MATRIX_1DOUBLE_OUT_MATRIX};
    char* test_function_in_1matrix_TYPE_out_matrix_name[] = {IN_1MATRIX_1DOUBLE_OUT_MATRIX_NAME};
    char* test_function_in_1matrix_TYPE_out_matrix_in0[] = {IN_1MATRIX_1DOUBLE_OUT_MATRIX_IN0};
    size = sizeof(test_function_in_1matrix_TYPE_out_matrix_in0)/sizeof(char*);
    matrix_t** test_function_in_1matrix_TYPE_out_matrix_in0_input = chartab2matrixtab(test_function_in_1matrix_TYPE_out_matrix_in0, size, data_path);
    TYPE test_function_in_1matrix_TYPE_out_matrix_in1[] = {IN_1MATRIX_1DOUBLE_OUT_MATRIX_IN1};
    char* test_function_in_1matrix_TYPE_out_matrix_out[] = {IN_1MATRIX_1DOUBLE_OUT_MATRIX_OUT};
    matrix_t** test_function_in_1matrix_TYPE_out_matrix_out_output = NULL;
    size = sizeof(test_function_in_1matrix_TYPE_out_matrix_out)/sizeof(char*);
    test_function_in_1matrix_TYPE_out_matrix_out_output = chartab2matrixtab(test_function_in_1matrix_TYPE_out_matrix_out, size, data_path);
    for (int i = 0; i < size; i++)
        if(!test_function_in_matrix_TYPE_out_matrix_f(test_function_in_1matrix_TYPE_out_matrix_f[i],test_function_in_1matrix_TYPE_out_matrix_in0_input[i], test_function_in_1matrix_TYPE_out_matrix_in1[i], test_function_in_1matrix_TYPE_out_matrix_out_output[i], test_function_in_1matrix_TYPE_out_matrix_name[i]))ret=1;
    free_matrixtab(test_function_in_1matrix_TYPE_out_matrix_in0_input, size);  
    free_matrixtab(test_function_in_1matrix_TYPE_out_matrix_out_output, size);   
    
    matrix_t*(*test_function_in_1matrix_1int_out_matrix_f[])(const matrix_t*, int) = {IN_1MATRIX_1INT_OUT_MATRIX};
    char* test_function_in_1matrix_1int_out_matrix_name[] = {IN_1MATRIX_1INT_OUT_MATRIX_NAME};
    char* test_function_in_1matrix_1int_out_matrix_in0[] = {IN_1MATRIX_1INT_OUT_MATRIX_IN0};
    size = sizeof(test_function_in_1matrix_1int_out_matrix_in0)/sizeof(char*);
    matrix_t** test_function_in_1matrix_1int_out_matrix_in0_input = chartab2matrixtab(test_function_in_1matrix_1int_out_matrix_in0, size, data_path);
    int test_function_in_1matrix_1int_out_matrix_in1[] = {IN_1MATRIX_1INT_OUT_MATRIX_IN1};
    char * test_function_in_1matrix_1int_out_matrix_out[] = {IN_1MATRIX_1INT_OUT_MATRIX_OUT};
    matrix_t** test_function_in_1matrix_1int_out_matrix_out_output = NULL;
    size = sizeof(test_function_in_1matrix_1int_out_matrix_out)/sizeof(char*);
    test_function_in_1matrix_1int_out_matrix_out_output = chartab2matrixtab(test_function_in_1matrix_1int_out_matrix_out, size, data_path);
    for (int i = 0; i < size; i++)
        if(!test_function_in_matrix_int_out_matrix_f(test_function_in_1matrix_1int_out_matrix_f[i],test_function_in_1matrix_1int_out_matrix_in0_input[i], test_function_in_1matrix_1int_out_matrix_in1[i], test_function_in_1matrix_1int_out_matrix_out_output[i], test_function_in_1matrix_1int_out_matrix_name[i]))ret=1; 
    free_matrixtab(test_function_in_1matrix_1int_out_matrix_in0_input, size);  
    free_matrixtab(test_function_in_1matrix_1int_out_matrix_out_output, size);  
    
    matrix_t*(*test_function_in_1matrix_2int_out_matrix_f[])(const matrix_t*, size_t, size_t) = {IN_1MATRIX_2INT_OUT_MATRIX};
    char* test_function_in_1matrix_2int_out_matrix_name[] = {IN_1MATRIX_2INT_OUT_MATRIX_NAME};
    char* test_function_in_1matrix_2int_out_matrix_in0[] = {IN_1MATRIX_2INT_OUT_MATRIX_IN0};
    size = sizeof(test_function_in_1matrix_2int_out_matrix_in0)/sizeof(char*);
    matrix_t** test_function_in_1matrix_2int_out_matrix_in0_input = chartab2matrixtab(test_function_in_1matrix_2int_out_matrix_in0, size, data_path);
    int test_function_in_1matrix_2int_out_matrix_in1[] = {IN_1MATRIX_2INT_OUT_MATRIX_IN1};
    int test_function_in_1matrix_2int_out_matrix_in2[] = {IN_1MATRIX_2INT_OUT_MATRIX_IN2};
    char * test_function_in_1matrix_2int_out_matrix_out[] = {IN_1MATRIX_2INT_OUT_MATRIX_OUT};
    matrix_t** test_function_in_1matrix_2int_out_matrix_out_output = NULL;
    size = sizeof(test_function_in_1matrix_2int_out_matrix_out)/sizeof(char*);
    test_function_in_1matrix_2int_out_matrix_out_output = chartab2matrixtab(test_function_in_1matrix_2int_out_matrix_out, size, data_path);
    for (int i = 0; i < size; i++)
        if(!test_function_in_matrix_2int_out_matrix_f(test_function_in_1matrix_2int_out_matrix_f[i],test_function_in_1matrix_2int_out_matrix_in0_input[i],test_function_in_1matrix_2int_out_matrix_in1[i], test_function_in_1matrix_2int_out_matrix_in2[i], test_function_in_1matrix_2int_out_matrix_out_output[i], test_function_in_1matrix_2int_out_matrix_name[i]))ret=1; 
    free_matrixtab(test_function_in_1matrix_2int_out_matrix_in0_input, size);  
    free_matrixtab(test_function_in_1matrix_2int_out_matrix_out_output, size); 
    
    if(!test_error_cases())ret=1;
    ret=0;
    return ret;
}