#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <sys/sysinfo.h>
#include "thread_pool.h"
#include "fifo.h"
thread_pool_t thread_pool;
char * format_time(const long long input_time, char* format)
{
    const char* formats[10] = {"a","d","h","m","s","ms","µs","ns","ps","fs"};
    const int  timescales[10] = {1,365,24,60,60,1000,1000,1000,1000,1000};
    const int width[10] = {0,3,2,2,2,3,3,3,3,3};
    long long timestamp[10];
    int  scale, i, j, k;
    size_t bufsz;
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

typedef struct{
    int work_len;
    int* ret;
}work_t;

int is_prime(int n)
{
    if(n==0 || n==1)return 0;
    for (int i=2; i<=floor(sqrt(n)); i++)
        if (n%i == 0)
            return 0;
    return 1;
}

void prime_calc(void *arg, int n)
{
    // Recherche de tous les nombres premiers entre work->min et work->max
    work_t *work = arg;
    int min = n*work->work_len;
    int max = min + work->work_len - 1;
    int j=0;
    for(int i = min; i<max; i++){
       if (is_prime(i))
            work->ret[min + j++] = i;
    }
}

static work_t _main_thread(int work_nb, int work_len)
{
    // Working without thread_pool
    printf("\t* le thread principal: ");
    fflush(stdout);
    work_t work = {work_len, calloc(work_nb*work_len, sizeof(int))};
    long long time = mstime();
    for (int i=0;i<work_nb;i++)
        prime_calc(&work, i);
    char *formatted_time = format_time(mstime() - time, "ms");
    printf("%s\n", formatted_time);
    free(formatted_time);
    return work;
}

static work_t _thread_pool(int work_nb, int work_len)
{
    // Working with thread_pool
    printf("\t* le pool de threads: ");
    fflush(stdout);
    work_t work = {work_len, calloc(work_nb*work_len, sizeof(int))};
    thread_pool_work_t thpool_work = {0, NULL, prime_calc, &work};
    long long time = mstime();
    for (int i=0;i<work_nb;i++){
        if(thread_pool_queue_work(&thread_pool, &thpool_work, i) != THREAD_POOL_OK){
            printf("\x1b[31mproblem1\x1b[0m\n");
        }
    }
    if(thread_pool_wait(&thread_pool) != THREAD_POOL_OK){
        printf("\x1b[31mproblem2\x1b[0m\n");
    }
    char *formatted_time = format_time(mstime() - time, "ms");
    printf("%s\n", formatted_time);
    free(formatted_time);
    return work;
}

int main(int argc, char ** argv)
{
    // "Traitement de 10 sacs de 2kg de données par 8 thread"
    int work_nb = 10;
    int work_len = 2000;
    int threads = 8;
    if (argc > 1)work_nb = atoi(argv[1]);
    if (argc > 2)work_len = atoi(argv[2]);
    if (argc > 3)threads = atoi(argv[3]);
    printf("This system has %d processors configured and %d processors available.\n", get_nprocs_conf(), get_nprocs());
        
    //Creating thread pool
    printf("Création d'un pool de \x1b[34m%d\x1b[0m threads\n", threads);
    if(thread_pool_create(&thread_pool, threads, NULL) != THREAD_POOL_OK){
        printf("\x1b[31mproblem0\x1b[0m\n");
        return 1;
    }
    
    //Preparing test functions
    int num_func = 2;
    work_t (*func_tab[2])(int work_nb, int work_len) = {_main_thread, _thread_pool};
    work_t res[num_func];
    //Executing test functions
    printf("Traitement de \x1b[34m%d\x1b[0m sacs de \x1b[34m%dg\x1b[0m de données par :\n", work_nb, work_len);
    for (int i=0; i<num_func; i++)
        res[i] = (func_tab[i])(work_nb, work_len);
    
    //Testing result
    printf("Vérification de l'intégrité des résultats: ");
    int ok =1;
    for (int i=1; i<num_func; i++){
        for (int j=0; j<work_nb; j++){
            for (int k=0; k<work_len; k++){
                if(res[i].ret[j*work_len+k] != res[0].ret[j*work_len+k]){
                    printf("\x1b[31mInconsistency in function %d\x1b[0m\n", i);
                    ok = 0;
                    break;
                }
            }
        }
    }     
    if(ok)printf("\x1b[32mOK\x1b[0m\n");
    // Freeing results
    for (int i=0; i<num_func; i++)
        free(res[i].ret);
    //Displaying the result
    // int len = snprintf(NULL,0, "[ ");
    // for (int i=0;i<work_nb;i++)
        // for(int j=0;j<work_len;j++)
            // if(works[i].ret[j] != 0)
                // len+= snprintf(NULL,0, "%d, ", works[i].ret[j]);
    // len += snprintf(NULL,0, " ]");
    // char *str = calloc(len, sizeof(*str));
    // char* ptr = str;
    // ptr += sprintf(ptr,"[ ");
    // for (int i=0;i<work_nb;i++)
        // for(int j=0;j<work_len;j++)
            // if(works[i].ret[j] != 0)
                // ptr += sprintf(ptr,"%d, ", works[i].ret[j]);
    // sprintf(ptr-2," ]");
    // printf("%s\n", str);
    // free(str);
    // for (int i=0;i<work_nb;i++)
        // free(works[i].ret);
    //Destroying thread pool
    if(thread_pool_destroy(&thread_pool) != THREAD_POOL_OK){
        printf("\x1b[31mproblem3\x1b[0m\n");
        return 1;
    }
    return 0;
}