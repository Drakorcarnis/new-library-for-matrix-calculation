#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include "includes/matrix.h"

#define DATA_PATH "/home/Drakorcarnis/new_library/src/data"
#define PRECISION 1e-15
#define INT_PRECISION 15
double precision = PRECISION;
double int_precision = INT_PRECISION;

#define IN_1MATRIX_OUT_DOUBLE       matrix_det_plu_f,   matrix_det_cholesky_f
#define IN_1MATRIX_OUT_DOUBLE_IN    "matrix",           "matrix_sym"
#define IN_1MATRIX_OUT_DOUBLE_OUT   1.50927e+09,        1
#define IN_1MATRIX_OUT_DOUBLE_NAME  "matrix_det_plu_f", "matrix_det_cholesky_f"

#define IN_1MATRIX_OUT_MATRIX       matrix_transp_f,    matrix_inverse_plu_f,   matrix_inverse_cholesky_f
#define IN_1MATRIX_OUT_MATRIX_IN    "matrix",           "matrix",               "matrix_sym"
#define IN_1MATRIX_OUT_MATRIX_OUT   "matrix_transp",    "matrix_inv",           "matrix_sym_inv"
#define IN_1MATRIX_OUT_MATRIX_NAME  "matrix_transp_f",  "matrix_inverse_plu_f", "matrix_inverse_cholesky_f"

#define IN_2MATRIX_OUT_MATRIX       matrix_add_f,   matrix_mult_f,      matrix_solve_plu_f,     matrix_solve_cholesky_f,  
#define IN_2MATRIX_OUT_MATRIX_IN0   "matrix",       "matrix",           "matrix",               "matrix_sym",                 
#define IN_2MATRIX_OUT_MATRIX_IN1   "matrix",       "matrix",           "B",                    "B",                      
#define IN_2MATRIX_OUT_MATRIX_OUT   "matrix_add",   "matrix_mult",      "X",                    "X_sym",                      
#define IN_2MATRIX_OUT_MATRIX_NAME  "matrix_add_f", "matrix_mult_f",    "matrix_solve_plu_f",   "matrix_sym_solve_cholesky_f",

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

#define IN_1MATRIX_2INT_OUT_MATRIX      matrix_shrink_f
#define IN_1MATRIX_2INT_OUT_MATRIX_IN0  "matrix"
#define IN_1MATRIX_2INT_OUT_MATRIX_IN1  1
#define IN_1MATRIX_2INT_OUT_MATRIX_IN2  2
#define IN_1MATRIX_2INT_OUT_MATRIX_OUT  "matrix_shrink"
#define IN_1MATRIX_2INT_OUT_MATRIX_NAME "matrix_shrink_f"

typedef struct {
    char *test_name;
    int result;
    long long utime;
}result_t;

char * format_time(const long long input_time, char* format)
{
    const char* formats[10] = {"a","d","h","m","s","ms","µs","ns","ps","fs"};
    const int  timescales[10] = {1,365,24,60,60,1000,1000,1000,1000,1000};
    const int width[10] = {0,3,2,2,2,3,3,3,3,3};
    long long timestamp[10];
    int  scale, i, j, count;
    char *ret = malloc(30*sizeof(char));
    for (scale = 1; (scale <= 10) && (strcmp(format, formats[scale-1]) != 0); scale++);
    if(scale > 10)return(NULL); // Unsupported format, you filthy rat !
    if(input_time<=0){sprintf(ret, "0%s\n", format);return(ret);} // Quickly handle case 0
    for (i = 0; i <= scale; i++)
    {
        timestamp[i] = input_time;
        for (j = scale - 1; j > i; j--)timestamp[i] = timestamp[i]/timescales[j];
        if(i > 0)timestamp[i] = timestamp[i]%timescales[i];
    }
    for (i=0; i < scale && timestamp[i] == 0; i++ );
    for (j=scale; j > i && timestamp[j] == 0; j-- );
    count = sprintf(ret, "%lld%s",timestamp[i],formats[i]);
    for (i=i+1; i < j && (count +=(int)sprintf(ret + count, "%0*lld%s",width[i],timestamp[i],formats[i])); i++ );
    return(ret);
}

long long int utime(void)
{
   struct timeval tv;
   gettimeofday(&tv, NULL);
   return((long long)(tv.tv_usec + (long long)tv.tv_sec * 1000000));
}

void process_result(result_t res)
{
    char *formatted_time = format_time(res.utime, "µs");
    if (res.result)
        printf("\033[0;32m[OK]");
    else
        printf("\033[0;31m[NOK]");
    printf("%s (%s)\033[0m\n", res.test_name, formatted_time);
    free(formatted_time);
}

int test_matrix_equality(const matrix_t *matrix1, const matrix_t *matrix2){
    // matrix_display(matrix1);printf("\n");matrix_display(matrix2);
    if((matrix1->rows != matrix2->rows) || (matrix1->columns != matrix2->columns))return 0;
    for (int i=0; i<matrix1->rows; i++){
        for (int j=0; j<matrix1->columns; j++){
            if(fabs(creal(matrix1->coeff[i][j] - matrix2->coeff[i][j])) > precision)return 0;
            if(fabs(cimag(matrix1->coeff[i][j] - matrix2->coeff[i][j])) > precision)return 0;
        }
    }
    return 1;
}

result_t test_function_in_matrix_out_double_f(double complex(*function)(const matrix_t*), matrix_t* matrix, double complex expected, char *test_name)
{
    long long time = utime();
    double complex ret = (function)(matrix);
    long long time2 = utime();
    result_t res = {test_name, fabs(-1+creal(ret)/creal(expected)) < 2e-3, time2 - time};
    process_result(res);
    if (!res.result) printf("expected: %f\ngot: %f\n",creal(expected), creal(ret));
    return(res);
}

result_t test_function_in_matrix_out_matrix_f(matrix_t*(*function)(const matrix_t*), matrix_t* matrix, matrix_t* expected, char *test_name)
{
    long long time = utime();
    matrix_t *ret = (function)(matrix);
    long long time2 = utime();
    result_t res = {test_name, test_matrix_equality(expected, ret), time2 - time};
    process_result(res);
    if (!res.result){printf("expected:\n");matrix_display_exact(expected, int_precision);printf("got:\n");matrix_display_exact(ret, int_precision);}
    matrix_free(ret);
    return(res);
}

result_t test_function_in_2matrix_out_matrix_f(matrix_t*(*function)(const matrix_t*, const matrix_t*), matrix_t* matrix1, matrix_t* matrix2, matrix_t* expected, char *test_name)
{
    long long time = utime();
    matrix_t *ret = (function)(matrix1, matrix2);
    long long time2 = utime();
    result_t res = {test_name, test_matrix_equality(expected, ret), time2 - time};
    process_result(res);
    if (!res.result){printf("expected:\n");matrix_display_exact(expected, int_precision);printf("got:\n");matrix_display_exact(ret, int_precision);}
    matrix_free(ret);
    return(res);
}

result_t test_function_in_matrix_double_out_matrix_f(matrix_t*(*function)(const matrix_t*, double complex), matrix_t* matrix, double complex val, matrix_t* expected, char *test_name)
{
    long long time = utime();
    matrix_t *ret = (function)(matrix, val);
    long long time2 = utime();
    result_t res = {test_name, test_matrix_equality(expected, ret), time2 - time};
    process_result(res);
if (!res.result){printf("expected:\n");matrix_display_exact(expected, int_precision);printf("got:\n");matrix_display_exact(ret, int_precision);}
    matrix_free(ret);
    return(res);
}

result_t test_function_in_matrix_int_out_matrix_f(matrix_t*(*function)(const matrix_t*, int), matrix_t* matrix, int val, matrix_t* expected, char *test_name)
{
    long long time = utime();
    matrix_t *ret = (function)(matrix, val);
    long long time2 = utime();
    result_t res = {test_name, test_matrix_equality(expected, ret), time2 - time};
    process_result(res);
if (!res.result){printf("expected:\n");matrix_display_exact(expected, int_precision);printf("got:\n");matrix_display_exact(ret, int_precision);}
    matrix_free(ret);
    return(res);
}

matrix_t** chartab2matrixtab(char ** filetab, ssize_t size, char *data_path)
{
    matrix_t** matrixtab = malloc(sizeof(matrix_t*) * size);
    for (int i = 0; i < size; i++)
    {
        ssize_t bufsz = snprintf(NULL, 0, "%s/%s.txt",data_path,filetab[i]);
        char* filename = malloc(bufsz + 1);
        snprintf(filename, bufsz + 1, "%s/%s.txt",data_path,filetab[i]);
        matrixtab[i] = file2matrix(filename);
        free(filename);
    }
    return matrixtab;
}

void free_matrixtab(matrix_t** matrixtab, ssize_t size)
{
    for (int i = 0; i < size; i++)free(matrixtab[i]);
    free(matrixtab);
    matrixtab = NULL;
}

int main(int argc, char **argv) {
    char *data_path = DATA_PATH;
    if(argc > 1){
        ssize_t bufsz = snprintf(NULL, 0, "1e-%s",argv[1]);
        char *str_precision = malloc(bufsz);
        snprintf(str_precision, bufsz + 1, "1e-%s",argv[1]);
        precision = strtod(str_precision, NULL);
        int_precision = atoi(argv[1]);
    }
    ssize_t size;
    
    double complex(*test_function_in_1matrix_out_double_f[])(const matrix_t*) = {IN_1MATRIX_OUT_DOUBLE};
    char* test_function_in_1matrix_out_double_name[] = {IN_1MATRIX_OUT_DOUBLE_NAME};
    char* test_function_in_1matrix_out_double_in[] = {IN_1MATRIX_OUT_DOUBLE_IN};
    size = sizeof(test_function_in_1matrix_out_double_in)/sizeof(char*);
    matrix_t** test_function_in_1matrix_out_double_in_input = chartab2matrixtab(test_function_in_1matrix_out_double_in, size, data_path);
    double complex test_function_in_1matrix_out_double_out[] = {IN_1MATRIX_OUT_DOUBLE_OUT}; 
    for (int i = 0; i < size; i++)
    test_function_in_matrix_out_double_f(test_function_in_1matrix_out_double_f[i],test_function_in_1matrix_out_double_in_input[i], test_function_in_1matrix_out_double_out[i], test_function_in_1matrix_out_double_name[i]);    
    free_matrixtab(test_function_in_1matrix_out_double_in_input, size);
    
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
    test_function_in_matrix_out_matrix_f(test_function_in_1matrix_out_matrix_f[i],test_function_in_1matrix_out_matrix_in_input[i], test_function_in_1matrix_out_matrix_out_output[i], test_function_in_1matrix_out_matrix_name[i]); 
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
    test_function_in_2matrix_out_matrix_f(test_function_in_2matrix_out_matrix[i],test_function_in_2matrix_out_matrix_in0_input[i], test_function_in_2matrix_out_matrix_in1_input[i], test_function_in_2matrix_out_matrix_out_output[i], test_function_in_2matrix_out_matrix_name[i]);  
    free_matrixtab(test_function_in_2matrix_out_matrix_in0_input, size);  
    free_matrixtab(test_function_in_2matrix_out_matrix_in1_input, size);  
    free_matrixtab(test_function_in_2matrix_out_matrix_out_output, size);  
    
    matrix_t*(*test_function_in_1matrix_double_out_matrix_f[])(const matrix_t*, double complex) = {IN_1MATRIX_1DOUBLE_OUT_MATRIX};
    char* test_function_in_1matrix_double_out_matrix_name[] = {IN_1MATRIX_1DOUBLE_OUT_MATRIX_NAME};
    char* test_function_in_1matrix_double_out_matrix_in0[] = {IN_1MATRIX_1DOUBLE_OUT_MATRIX_IN0};
    size = sizeof(test_function_in_1matrix_double_out_matrix_in0)/sizeof(char*);
    matrix_t** test_function_in_1matrix_double_out_matrix_in0_input = chartab2matrixtab(test_function_in_1matrix_double_out_matrix_in0, size, data_path);
    double test_function_in_1matrix_double_out_matrix_in1[] = {IN_1MATRIX_1DOUBLE_OUT_MATRIX_IN1};
    char* test_function_in_1matrix_double_out_matrix_out[] = {IN_1MATRIX_1DOUBLE_OUT_MATRIX_OUT};
    matrix_t** test_function_in_1matrix_double_out_matrix_out_output = NULL;
    size = sizeof(test_function_in_1matrix_double_out_matrix_out)/sizeof(char*);
    test_function_in_1matrix_double_out_matrix_out_output = chartab2matrixtab(test_function_in_1matrix_double_out_matrix_out, size, data_path);
    for (int i = 0; i < size; i++)
    test_function_in_matrix_double_out_matrix_f(test_function_in_1matrix_double_out_matrix_f[i],test_function_in_1matrix_double_out_matrix_in0_input[i], test_function_in_1matrix_double_out_matrix_in1[i], test_function_in_1matrix_double_out_matrix_out_output[i], test_function_in_1matrix_double_out_matrix_name[i]);
    free_matrixtab(test_function_in_1matrix_double_out_matrix_in0_input, size);  
    free_matrixtab(test_function_in_1matrix_double_out_matrix_out_output, size);   
    
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
    test_function_in_matrix_int_out_matrix_f(test_function_in_1matrix_1int_out_matrix_f[i],test_function_in_1matrix_1int_out_matrix_in0_input[i], test_function_in_1matrix_1int_out_matrix_in1[i], test_function_in_1matrix_1int_out_matrix_out_output[i], test_function_in_1matrix_1int_out_matrix_name[i]); 
    free_matrixtab(test_function_in_1matrix_1int_out_matrix_in0_input, size);  
    free_matrixtab(test_function_in_1matrix_1int_out_matrix_out_output, size);  
    
    return 0;
}