#define _XOPEN_SOURCE 700
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <sys/time.h>
#include "matrix.h"
#include "matrix_tools.h"

static int sanity_check(const void *pointer, const char *function_name)
{
    if(!pointer){
        fprintf(stderr, "%s: NULL pointer\n",function_name);
        return 0;
    }
    return 1;
}

static double str2double(const char *str, unsigned int len)
{
    char *num = calloc(len+1, sizeof(char));
    char *den = calloc(len+1, sizeof(char));
    if ((!num) || (!den)){
        perror(__func__);
        fprintf(stderr, "Wanna alloc %d\n", len+1); 
        return -1;
    }
    den[0] = '1';
    for (unsigned int i = 0; i < len; i++) {
        if (str[i] == '/'){
            memcpy(den, str+i+1, len-i-1);
            break;
        } else {
            num[i] = str[i];
        }
    }
    double ret = strtod(num,NULL)/strtod(den,NULL);
    free(num);
    free(den);
    return ret;
}

static void _matrix_display(const matrix_t *matrix, int precision, FILE *stream)
{
    double val;
    int len = 0;
    for (unsigned int i = 0; i < matrix->rows; i++) {
        if (!precision)
            fprintf(stream, "[ ");
        for (unsigned int j = 0; j < matrix->columns; j++){
            val = matrix->coeff[i][j];
            if (precision){
                fprintf(stream, "%.*g ", precision, val);
                continue;
            }
            if(fabs(val) < 1e-10) {
                len = snprintf(NULL, 0, "0");
                fprintf(stream, "0");
            } else {
                len = snprintf(NULL, 0, "%.3g", val);
                fprintf(stream, "%.3g", val);
            }
            for (int k=0; k< 10-len; k++)fprintf(stream, " ");
        }
        if (!precision)
            fprintf(stream, "]");
        fprintf(stream, "\n");
    }
}

matrix_t * matrix_random(int rows, int columns)
{
    matrix_t *matrix = matrix_create(rows, columns);
    srand((unsigned int)time(NULL));
    for (unsigned int i = 0; i < matrix->rows; i++) {
        for (unsigned int j = 0; j < matrix->columns; j++) {
            matrix->coeff[i][j] = (double)(pow(-1.0, rand())*(rand()%10));
        }
    }
    return matrix;
}

matrix_t * matrix_symetric_random(int rows, int columns)
{
    double random;
    matrix_t *matrix = matrix_create(rows, columns);
    srand((unsigned int)time(NULL));
    for (unsigned int i = 0; i < matrix->rows; i++) {
        for (unsigned int j = 0; j <= i; j++) {
            random = (double)(pow(-1.0, rand())*(rand()%10));
            matrix->coeff[i][j] = matrix->coeff[j][i] = random;
        }
    }
    return matrix;
}

matrix_t * str2matrix(int argc, char **argv, char separator)
{
    unsigned int i, j, columns = 0, len, trigger, flag = 0, count, size;
    matrix_t *matrix = NULL;
    for (i = 0; i < (unsigned int)argc; i++) {
        count = 0;
        for (j = 0; argv[i][j]; j++){
            if(argv[i][j] != separator && argv[i][j] != '\t' && argv[i][j] != '\n' && argv[i][j] != '\r'){
                if(!flag)
                    count++;
                flag = 1;
            } else {
                flag = 0;
            }
        }
        if (count>columns)columns=count;
    }
    matrix = matrix_create(argc, columns);
    for (i = 0; i < matrix->rows; i++) {
        size = strlen(argv[i]);
        len = 0;
        flag = 0;
        trigger = 0;
            for (j = 0; j <= size; j++) {
                if(argv[i][j] != separator && argv[i][j] != '\t' && argv[i][j] != '\n' && argv[i][j] != '\r'){
                    if(!flag)trigger=j;
                    flag = 1;
                } else {
                    if(flag)matrix->coeff[i][len++]=str2double(argv[i]+trigger,j-trigger);
                    flag = 0;
                }
            }
        if(trigger <= size && flag && len < columns)matrix->coeff[i][len]=str2double(argv[i]+trigger,size-trigger);
    }
    return matrix;
}

int matrix2file(matrix_t *matrix, char * filename)
{
    FILE *fp = fopen(filename, "w");
    if(!fp){
        perror(__func__);
        return 0;
    }
    _matrix_display(matrix, 15, fp);
    fclose(fp);
    return 1;
}
    
matrix_t * file2matrix(char *filename)
{
    FILE *fp = fopen(filename, "r");
    if(!fp){
        perror(__func__);
        return NULL;
    }
    int argc = 0, i;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    matrix_t *matrix = NULL;
    char **argv = malloc(sizeof(char*));
    if(!argv){
        perror(__func__);
        return NULL;
    }
    while ((read = getline(&line, &len, fp)) != -1) {
        if (read > 1 || ((read == 1) && (strcmp(line, "\n") != 0))){
            argv = realloc(argv, (argc + 1) * sizeof(char*));
            if(!argv){
                perror(__func__);
                goto finally;
            }
            argv[argc] = calloc(read+1, sizeof(char));
            if(!argv[argc]){
                perror(__func__);
                goto finally;
            }
            strncpy(argv[argc++], line, read);
        }    
    }
    matrix = str2matrix(argc, argv, ' ');
finally:
    for (i=0; i<argc; i++)
    {
        free(argv[i]);
    }
    free(argv);
    if(line)free(line);
    fclose(fp);
    return matrix;
}

// Matrix display functions

void matrix_display(const matrix_t *matrix)
{
    _matrix_display(matrix, 0, stdout);
}
void matrix_display_exact(const matrix_t *matrix, int precision)
{
    _matrix_display(matrix, precision, stdout);
}

char * format_time(const long long input_time, char* format)
{
    const char* formats[10] = {"a","d","h","m","s","ms","µs","ns","ps","fs"};
    const int  timescales[10] = {1,365,24,60,60,1000,1000,1000,1000,1000};
    const int width[10] = {0,3,2,2,2,3,3,3,3,3};
    long long timestamp[10];
    int  scale, i, j, k;
    ssize_t bufsz;
    char *ret = NULL;
    for (scale = 1; (scale <= 10) && (strcmp(format, formats[scale-1]) != 0); scale++);
    if(scale > 10)return(NULL); // Unsupported format, you filthy rat !
    if(input_time<=0){ // Quickly handle case 0
        bufsz = snprintf(NULL, 0, "0%s", format);
        ret = malloc((bufsz+1)*sizeof(*ret));
        if(!ret){perror(__func__);}
        snprintf(ret, bufsz, "0%s", format);
        return(ret);
    }
    for (i = 0; i <= scale; i++){
        timestamp[i] = input_time;
        for (j = scale - 1; j > i; j--)timestamp[i] = timestamp[i]/timescales[j];
        if(i > 0)timestamp[i] = timestamp[i]%timescales[i];
    }
    for (i=0; i < scale && timestamp[i] == 0; i++);
    for (j=scale; j > i && timestamp[j] == 0; j--);
    bufsz = snprintf(NULL, 0, "%lld%s",timestamp[i],formats[i]);
    for (k=i+1; k < j && (bufsz += snprintf(NULL, 0, "%0*lld%s",width[k],timestamp[k],formats[k])); k++ );
    ret = malloc((bufsz+1)*sizeof(*ret));
    if(!ret){perror(__func__);}
    bufsz = sprintf(ret, "%lld%s",timestamp[i],formats[i]);
    for (k=i+1; k < j && (bufsz +=(int)sprintf(ret + bufsz, "%0*lld%s",width[k],timestamp[k],formats[k])); k++ );
    return(ret);
}

long long mstime(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return((long long)(tv.tv_usec/1000 + (long long)tv.tv_sec * 1000));
    // return((long long)(1e3*clock()/CLOCKS_PER_SEC));
}

int test_matrix_equality(const matrix_t *matrix1, const matrix_t *matrix2, int precision)
{
    if(!sanity_check(matrix1, __func__))return 0; 
    if(!sanity_check(matrix2, __func__))return 0; 
    if((matrix1->rows != matrix2->rows) || (matrix1->columns != matrix2->columns))return 0;
    for (unsigned int i=0; i<matrix1->rows; i++){
        for (unsigned int j=0; j<matrix1->columns; j++){
            if(fabs(matrix1->coeff[i][j] != matrix1->coeff[i][j]) )return 0;// These checks are necessary 
            if(fabs(matrix2->coeff[i][j] != matrix2->coeff[i][j]))return 0;//  in case of nan values
            if(fabs(matrix1->coeff[i][j] - matrix2->coeff[i][j]) > pow(10, -precision))return 0;
        }
    }
    return 1;
}