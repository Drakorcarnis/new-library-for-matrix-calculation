#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include "includes/matrix.h"


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

int test_matrix_add_f(const matrix_t *matrix1, const matrix_t *matrix2)
{
    long long time = utime();
    matrix_t *sum_matrix = matrix_add_f(matrix1, matrix2);
    long long time2 = utime();
    if(!sum_matrix)return(0);
    char *formatted_time = format_time(time2 - time, "µs");
    printf("\n(%s)M+M =\n", formatted_time);
    matrix_display(sum_matrix);
    matrix_free(sum_matrix);
    free(formatted_time);
    return(1);
}

int test_matrix_mult_scalar_f(const matrix_t *matrix, double complex lambda)
{
    long long time = utime();
    matrix_t *mult_matrix = matrix_mult_scalar_f(matrix, lambda);
    long long time2 = utime();
    if(!mult_matrix)return(0);
    char *formatted_time = format_time(time2 - time, "µs");
    printf("\n(%s)5*M =\n", formatted_time);
    matrix_display(mult_matrix);
    matrix_free(mult_matrix);
    free(formatted_time);
    return(1);
}

int test_matrix_mult_f(const matrix_t *matrix1, const matrix_t *matrix2)
{
    long long time = utime();
    matrix_t *mat_mult_matrix = matrix_mult_f(matrix1, matrix2);
    long long time2 = utime();
    if(!mat_mult_matrix)return(0);
    char *formatted_time = format_time(time2 - time, "µs");
    printf("\n(%s)M1*M2 =\n", formatted_time);
    matrix_display(mat_mult_matrix);
    matrix_free(mat_mult_matrix);
    free(formatted_time);
    return(1);
}

int test_matrix_pow_f(const matrix_t *matrix, int pow)
{
    long long time = utime();
    matrix_t *pow_matrix = matrix_pow_f(matrix, pow);
    long long time2 = utime();
    if(!pow_matrix)return(0);
    char *formatted_time = format_time(time2 - time, "µs");
    printf("\n(%s)M^3 =\n", formatted_time);
    matrix_display(pow_matrix);
    matrix_free(pow_matrix);
    free(formatted_time);
    return(1);
}

int test_matrix_solve_raw_f(const matrix_t *matrix1, matrix_t *matrix2)
{
    long long time = utime();
    matrix_t *matrix_cramer = matrix_solve_raw_f(matrix1, matrix2);
    long long time2 = utime();
    if(!matrix_cramer)return(0);
    char *formatted_time = format_time(time2 - time, "µs");
    printf("\n(%s)raw X =\n", formatted_time);
    matrix_display(matrix_cramer);
    matrix_free(matrix_cramer);
    free(formatted_time);
    return(1);
}

int test_matrix_solve_plu_f(matrix_t *matrix1, matrix_t *matrix2)
{
    long long time = utime();
    matrix_t *matrix_cramer = matrix_solve_plu_f(matrix1, matrix2);
    long long time2 = utime();
    if(!matrix_cramer)return(0);
    char *formatted_time = format_time(time2 - time, "µs");
    printf("\n(%s)plu X =\n", formatted_time);
    matrix_display(matrix_cramer);
    matrix_free(matrix_cramer);
    free(formatted_time);
    return(1);
}

int test_matrix_solve_cholesky_f(matrix_t *matrix1, matrix_t *matrix2)
{
    long long time = utime();
    matrix_t *matrix_cramer = matrix_solve_cholesky_f(matrix1, matrix2);
    long long time2 = utime();
    if(!matrix_cramer)return(0);
    char *formatted_time = format_time(time2 - time, "µs");
    printf("\n(%s)cho X =\n", formatted_time);
    matrix_display(matrix_cramer);
    matrix_free(matrix_cramer);
    free(formatted_time);
    return(1);
}

int test_matrix_det_raw_f(matrix_t *matrix)
{
    long long time = utime();
    double complex det = matrix_det_raw_f(matrix);
    long long time2 = utime();
    char *formatted_time = format_time(time2 - time, "µs");
    printf("\n(%s)raw |M| = %g\n", formatted_time, creal(det));
    return(1);
}
int test_matrix_det_plu_f(matrix_t *matrix)
{
    long long time = utime();
    double complex det = matrix_det_plu_f(matrix);
    long long time2 = utime();
    char *formatted_time = format_time(time2 - time, "µs");
    printf("\n(%s)plu |M| = %g\n", formatted_time, creal(det));
    free(formatted_time);
    return(1);
}
int test_matrix_det_cholesky_f(matrix_t *matrix)
{
    long long time = utime();
    double complex det = matrix_det_cholesky_f(matrix);
    long long time2 = utime();
    char *formatted_time = format_time(time2 - time, "µs");
    printf("\n(%s)cho |M| = %g\n", formatted_time, creal(det));
    free(formatted_time);
    return(1);
}
    
int test_function_f(matrix_t*(*function)(const matrix_t*), matrix_t* matrix, char* txt)
{
    long long time = utime();
    matrix_t *ret = (*function)(matrix);
    long long time2 = utime();
    char *formatted_time = format_time(time2 - time, "µs");
    if(!matrix)return(0);
    printf("\n(%s)%s =\n", formatted_time, txt);
    matrix_display(ret);
    matrix_free(ret);
    free(formatted_time);
    return(1);
}

int main(int argc, char *argv[]) {
    // matrix_t *matrix_p = str2matrix(argc-1, &argv[1], ',');
    int rank;
    if (argc < 2 || atoi(argv[1]) < 1)
        rank = 10;
    else
        rank = atoi(argv[1]);
    matrix_t *matrix1 = matrix_random(rank,rank);
    matrix_t *matrix2 = matrix_symetric_random(rank,rank);
    matrix_t *matrix3 = matrix_random(rank,1);
    if(!matrix1)return(-1);
    if(!matrix2)return(-1);
    if(!matrix3)return(-1);
    printf("M1 =\n");
    matrix_display(matrix1);
    printf("M2 =\n");
    matrix_display(matrix2);
    printf("M3 =\n");
    matrix_display(matrix3);
    
    test_function_f(&matrix_transp_f,matrix1, "transp(M)");
    test_matrix_add_f(matrix1, matrix2);
    test_matrix_mult_scalar_f(matrix1, 0.5);
    test_matrix_pow_f(matrix1, 3);
    test_matrix_mult_f(matrix1, matrix2);
    
    // test_function_f(&matrix_det_raw_f,matrix1, "|M|");
    // test_function_f(&matrix_com_f,matrix1, "com(M)");
    // test_function_f(&matrix_comp_f,matrix1, "comp(M)");
    // test_function_f(&matrix_inverse_raw_f, matrix1, "1/M");
    // test_matrix_det_raw_f(matrix2);
    
    test_matrix_solve_plu_f(matrix1, matrix3);
    test_function_f(&matrix_inverse_plu_f,matrix1,"plu 1/M");
    test_matrix_solve_cholesky_f(matrix2, matrix3);
    test_function_f(&matrix_inverse_plu_f,matrix2,"plu 1/M2");
    test_function_f(&matrix_inverse_cholesky_f,matrix2,"cho 1/M2");
    test_matrix_det_plu_f(matrix2);
    test_matrix_det_cholesky_f(matrix2);
    matrix_free(matrix1);
    matrix_free(matrix2);
    matrix_free(matrix3);
    return 0;
}