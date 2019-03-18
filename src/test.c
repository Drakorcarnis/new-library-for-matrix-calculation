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

int test_matrix_mult_scalar_f(const matrix_t *matrix, double lambda)
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
    printf("\n(%s)X =\n", formatted_time);
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
    printf("\n(%s)X =\n", formatted_time);
    matrix_display(matrix_cramer);
    matrix_free(matrix_cramer);
    free(formatted_time);
    return(1);
}

int test_matrix_plu_f(matrix_t *matrix)
{
    long long time = utime();
    plu_t *plu = matrix_plu_f(matrix);
    long long time2 = utime();
    char *formatted_time = format_time(time2 - time, "µs");
    if(!plu)return(0);
    printf("\nPLU computed in %s\n", formatted_time);
    printf("\nA =\n");
    matrix_display(matrix);
    printf("\nP =\n");
    matrix_display(plu->P);
    printf("\nL =\n");
    matrix_display(plu->L);
    printf("\nU =\n");
    matrix_display(plu->U);
    printf("\nP*L*U =\n");
    matrix_t *LU = matrix_mult_f(plu->L, plu->U);
    matrix_t *PLU = matrix_mult_f(plu->P, LU);
    matrix_display(PLU);
    matrix_free(LU);
    matrix_free(PLU);
    plu_free(plu);
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

int main() {
    // matrix_t *matrix_p = str2matrix(argc-1, &argv[1], ',');
    matrix_t *matrix1 = file2matrix("/home/Drakorcarnis/new_library/src/matrix1.txt");
    matrix_t *matrix2 = file2matrix("/home/Drakorcarnis/new_library/src/matrix2.txt");
    matrix_t *matrix3 = file2matrix("/home/Drakorcarnis/new_library/src/matrix3.txt");
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
    test_matrix_pow_f(matrix1, 10);
    // test_function_f(&matrix_det_raw_f,matrix1, "|M|");
    // test_function_f(&matrix_com_f,matrix1, "com(M)");
    // test_function_f(&matrix_comp_f,matrix1, "comp(M)");
    // test_function_f(&matrix_inverse_raw_f, matrix1, "1/M");
    test_matrix_mult_f(matrix1, matrix2);
    test_matrix_solve_plu_f(matrix1, matrix3);
    test_matrix_plu_f(matrix1);
    // test_function_f(&matrix_det_plu_f,matrix1, "|M|");
    test_function_f(&matrix_inverse_plu_f,matrix1,"1/M");
    matrix_free(matrix1);
    matrix_free(matrix2);
    matrix_free(matrix3);
    return 0;
}