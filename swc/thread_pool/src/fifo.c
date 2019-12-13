#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <signal.h>
#include <unistd.h>
#include <errno.h>
#include "fifo.h"

int fifo_init(fifo_t *fifo, size_t size)
{
    if(!fifo)return FIFO_FAIL_UNALLOCATED;
    int ret = FIFO_FAIL_MUTEX;
    fifo->size = size;
    fifo->buf = calloc(size, sizeof(void *));
    fifo->index_buf = malloc(size*sizeof(int));
    if(!fifo->buf){
        goto fail_alloc;
    }
    if(pthread_mutex_init(&fifo->curr_nb_mutex, NULL) != 0)goto fail_nb_mutex;
    if(pthread_mutex_init(&fifo->pop_mutex, NULL) != 0)goto fail_pop_mutex;
    if(pthread_mutex_init(&fifo->push_cond_mutex, NULL) != 0)goto fail_push_cond_mutex;
    if(pthread_mutex_init(&fifo->pop_cond_mutex, NULL) != 0)goto fail_pop_cond_mutex;
    if(pthread_mutex_init(&fifo->empty_cond_mutex, NULL) != 0)goto fail_empty_cond_mutex;
    if(pthread_cond_init(&fifo->push_cond, NULL) != 0)goto fail_push_cond;
    if(pthread_cond_init(&fifo->pop_cond, NULL) != 0)goto fail_pop_cond;
    if(pthread_cond_init(&fifo->empty_cond, NULL) != 0)goto fail_empty_cond;
    fifo->push_index = 0;
    fifo->pop_index = 0;
    fifo->curr_nb_elt = 0;
    goto success;
fail_empty_cond:
    ret = FIFO_FAIL_COND;
    pthread_cond_destroy(&fifo->pop_cond);
fail_pop_cond:
    ret = FIFO_FAIL_COND;
    pthread_cond_destroy(&fifo->push_cond);
fail_push_cond:
    pthread_mutex_destroy(&fifo->empty_cond_mutex);
fail_empty_cond_mutex:
    pthread_mutex_destroy(&fifo->pop_cond_mutex);
fail_pop_cond_mutex:
    pthread_mutex_destroy(&fifo->push_cond_mutex);
fail_push_cond_mutex:
    pthread_mutex_destroy(&fifo->pop_mutex);
fail_pop_mutex:
    pthread_mutex_destroy(&fifo->curr_nb_mutex);
fail_nb_mutex:
    free(fifo->buf);
    return ret;
fail_alloc:
    return FIFO_FAIL_MEM;
success:
    return FIFO_SUCCESS;
}

int fifo_push(fifo_t *fifo, void *elt, int wait)
{
    if(!fifo)return FIFO_FAIL_UNALLOCATED;
    if(pthread_mutex_lock(&fifo->push_cond_mutex))return FIFO_FAIL_MUTEX;
    if(fifo->curr_nb_elt == fifo->size && wait == FIFO_NO_WAIT){
        if(pthread_mutex_unlock(&fifo->push_cond_mutex)!=0)return FIFO_FAIL_MUTEX;
        return FIFO_FULL;
    }
    if(fifo->curr_nb_elt == fifo->size && wait == FIFO_WAIT){
        // printf("\x1b[31mPUSH WAITING (curr_nb_elt= %d, push_index = %d)\x1b[0m\n", fifo->curr_nb_elt, fifo->push_index);
        pthread_cond_wait(&fifo->push_cond, &fifo->push_cond_mutex);
    }
    fifo->buf[fifo->push_index] = elt;
    fifo->push_index = fifo->push_index+1 < fifo->size ? fifo->push_index+1:0;
    if(pthread_mutex_unlock(&fifo->push_cond_mutex))return FIFO_FAIL_MUTEX;
    if(pthread_mutex_lock(&fifo->pop_cond_mutex))return FIFO_FAIL_MUTEX;
    if(pthread_mutex_lock(&fifo->curr_nb_mutex))return FIFO_FAIL_MUTEX;
    fifo->curr_nb_elt++;
    if(pthread_mutex_unlock(&fifo->curr_nb_mutex))return FIFO_FAIL_MUTEX;
    if(pthread_cond_signal(&fifo->pop_cond))return FIFO_FAIL_COND;
    // printf("PUSHED work %p at index %d, curr_nb_elt = %d\n", elt, push_index, fifo->curr_nb_elt);
    if(pthread_mutex_unlock(&fifo->pop_cond_mutex))return FIFO_FAIL_MUTEX;
    return FIFO_SUCCESS;
}
int fifo_push_index(fifo_t *fifo, void *elt, int index, int wait)
{
    if(!fifo)return FIFO_FAIL_UNALLOCATED;
    if(pthread_mutex_lock(&fifo->push_cond_mutex))return FIFO_FAIL_MUTEX;
    if(fifo->curr_nb_elt == fifo->size && wait == FIFO_NO_WAIT){
        if(pthread_mutex_unlock(&fifo->push_cond_mutex)!=0)return FIFO_FAIL_MUTEX;
        return FIFO_FULL;
    }
    if(fifo->curr_nb_elt == fifo->size && wait == FIFO_WAIT){
        // printf("\x1b[31mPUSH WAITING (curr_nb_elt= %d, push_index = %d)\x1b[0m\n", fifo->curr_nb_elt, fifo->push_index);
        pthread_cond_wait(&fifo->push_cond, &fifo->push_cond_mutex);
    }
    fifo->buf[fifo->push_index] = elt;
    fifo->index_buf[fifo->push_index] = index;
    fifo->push_index = fifo->push_index+1 < fifo->size ? fifo->push_index+1:0;
    if(pthread_mutex_unlock(&fifo->push_cond_mutex))return FIFO_FAIL_MUTEX;
    if(pthread_mutex_lock(&fifo->pop_cond_mutex))return FIFO_FAIL_MUTEX;
    if(pthread_mutex_lock(&fifo->curr_nb_mutex))return FIFO_FAIL_MUTEX;
    fifo->curr_nb_elt++;
    if(pthread_mutex_unlock(&fifo->curr_nb_mutex))return FIFO_FAIL_MUTEX;
    if(pthread_cond_signal(&fifo->pop_cond))return FIFO_FAIL_COND;
    // printf("PUSHED work %p at index %d, curr_nb_elt = %d\n", elt, push_index, fifo->curr_nb_elt);
    if(pthread_mutex_unlock(&fifo->pop_cond_mutex))return FIFO_FAIL_MUTEX;
    return FIFO_SUCCESS;
}

int fifo_pop(fifo_t *fifo, void **elt, int wait)
{
    if(!fifo)return FIFO_FAIL_UNALLOCATED;
    if(pthread_mutex_lock(&fifo->pop_mutex)!=0)return FIFO_FAIL_MUTEX;
    if(wait == FIFO_NO_WAIT && fifo->curr_nb_elt == 0){
        if(pthread_mutex_unlock(&fifo->pop_mutex)!=0)return FIFO_FAIL_MUTEX;
        return FIFO_EMPTY;
    }
    if(pthread_mutex_lock(&fifo->pop_cond_mutex))return FIFO_FAIL_MUTEX;
    if(wait == FIFO_WAIT && fifo->curr_nb_elt == 0){
        printf("\x1b[31mPOP WAITING\x1b[0m\n");
        pthread_cond_wait(&fifo->pop_cond, &fifo->pop_cond_mutex);
    }
    if(pthread_mutex_unlock(&fifo->pop_cond_mutex))return FIFO_FAIL_MUTEX;
    *elt = fifo->buf[fifo->pop_index];
    fifo->buf[fifo->pop_index] = NULL;
    if(pthread_mutex_lock(&fifo->push_cond_mutex))return FIFO_FAIL_MUTEX;
    if(pthread_cond_signal(&fifo->push_cond))return FIFO_FAIL_COND;
    if(pthread_mutex_lock(&fifo->curr_nb_mutex))return FIFO_FAIL_MUTEX;
    fifo->curr_nb_elt--;
    if(pthread_mutex_unlock(&fifo->curr_nb_mutex))return FIFO_FAIL_MUTEX;
    // printf("POPPED work %p at index %d\n", *elt, fifo->pop_index);
    if(pthread_mutex_unlock(&fifo->push_cond_mutex))return FIFO_FAIL_MUTEX;
    fifo->pop_index = fifo->pop_index+1 < fifo->size ? fifo->pop_index+1:0;
    if(fifo->curr_nb_elt == 0){
        pthread_mutex_lock(&fifo->empty_cond_mutex);
        if(pthread_cond_signal(&fifo->empty_cond)!=0)return FIFO_FAIL_COND;
        pthread_mutex_unlock(&fifo->empty_cond_mutex);
    }
    if(pthread_mutex_unlock(&fifo->pop_mutex))return FIFO_FAIL_MUTEX;
    return FIFO_SUCCESS;
}
int fifo_pop_index(fifo_t *fifo, void **elt, int *index, int wait)
{
    if(!fifo)return FIFO_FAIL_UNALLOCATED;
    if(pthread_mutex_lock(&fifo->pop_mutex)!=0)return FIFO_FAIL_MUTEX;
    if(wait == FIFO_NO_WAIT && fifo->curr_nb_elt == 0){
        if(pthread_mutex_unlock(&fifo->pop_mutex)!=0)return FIFO_FAIL_MUTEX;
        return FIFO_EMPTY;
    }
    if(pthread_mutex_lock(&fifo->pop_cond_mutex))return FIFO_FAIL_MUTEX;
    if(wait == FIFO_WAIT && fifo->curr_nb_elt == 0){
        printf("\x1b[31mPOP WAITING\x1b[0m\n");
        pthread_cond_wait(&fifo->pop_cond, &fifo->pop_cond_mutex);
    }
    if(pthread_mutex_unlock(&fifo->pop_cond_mutex))return FIFO_FAIL_MUTEX;
    *elt = fifo->buf[fifo->pop_index];
    *index = fifo->index_buf[fifo->pop_index];
    fifo->buf[fifo->pop_index] = NULL;
    if(pthread_mutex_lock(&fifo->push_cond_mutex))return FIFO_FAIL_MUTEX;
    if(pthread_cond_signal(&fifo->push_cond))return FIFO_FAIL_COND;
    if(pthread_mutex_lock(&fifo->curr_nb_mutex))return FIFO_FAIL_MUTEX;
    fifo->curr_nb_elt--;
    if(pthread_mutex_unlock(&fifo->curr_nb_mutex))return FIFO_FAIL_MUTEX;
    // printf("POPPED work %p at index %d\n", *elt, fifo->pop_index);
    if(pthread_mutex_unlock(&fifo->push_cond_mutex))return FIFO_FAIL_MUTEX;
    fifo->pop_index = fifo->pop_index+1 < fifo->size ? fifo->pop_index+1:0;
    if(fifo->curr_nb_elt == 0){
        pthread_mutex_lock(&fifo->empty_cond_mutex);
        if(pthread_cond_signal(&fifo->empty_cond)!=0)return FIFO_FAIL_COND;
        pthread_mutex_unlock(&fifo->empty_cond_mutex);
    }
    if(pthread_mutex_unlock(&fifo->pop_mutex))return FIFO_FAIL_MUTEX;
    return FIFO_SUCCESS;
}

int fifo_wait_empty(fifo_t *fifo)
{
    if(!fifo)return FIFO_FAIL_UNALLOCATED;
    if(pthread_mutex_lock(&fifo->empty_cond_mutex))return FIFO_FAIL_MUTEX;
    if(fifo->curr_nb_elt != 0 && pthread_cond_wait(&fifo->empty_cond, &fifo->empty_cond_mutex)!=0){
        pthread_mutex_unlock(&fifo->empty_cond_mutex);
        return FIFO_FAIL_COND;
    }
    if(pthread_mutex_unlock(&fifo->empty_cond_mutex))return FIFO_FAIL_MUTEX;
    return FIFO_SUCCESS;
}
int fifo_current_size(fifo_t *fifo, int* curr_nb_elt)
{
    if(!fifo)return FIFO_FAIL_UNALLOCATED;
    *curr_nb_elt = fifo->curr_nb_elt;
    return FIFO_SUCCESS;
}
int fifo_destroy(fifo_t *fifo)
{
    if(!fifo)return FIFO_FAIL_UNALLOCATED;
    int ret = FIFO_SUCCESS;
    if(pthread_mutex_destroy(&fifo->curr_nb_mutex) != 0)ret = FIFO_FAIL_MUTEX;
    if(pthread_mutex_destroy(&fifo->pop_mutex) != 0)ret = FIFO_FAIL_MUTEX;
    if(pthread_mutex_destroy(&fifo->push_cond_mutex) != 0)ret = FIFO_FAIL_MUTEX;
    if(pthread_mutex_destroy(&fifo->pop_cond_mutex) != 0)ret = FIFO_FAIL_MUTEX;
    if(pthread_mutex_destroy(&fifo->empty_cond_mutex) != 0)ret = FIFO_FAIL_MUTEX;
    if(pthread_cond_destroy(&fifo->push_cond) != 0)ret = FIFO_FAIL_COND;
    if(pthread_cond_destroy(&fifo->pop_cond) != 0)ret = FIFO_FAIL_COND;
    if(pthread_cond_destroy(&fifo->empty_cond) != 0)ret = FIFO_FAIL_COND;
    free(fifo->buf);
    return ret;
}
int fifo_show(fifo_t *fifo)
{
    if(!fifo)return FIFO_FAIL_UNALLOCATED;
    char *string = NULL;
    size_t new_size, size = 0;
    printf("\x1b[34mNB OF ELEMENTS: %d\x1b[0m\n", fifo->curr_nb_elt);
    printf("\x1b[34mPUSH INDEX: %d\x1b[0m\n", fifo->push_index);
    printf("\x1b[34mPOP INDEX: %d\x1b[0m\n", fifo->pop_index);
    if(fifo->curr_nb_elt == 0)return FIFO_SUCCESS;
    for(size_t i=0; i<fifo->size; i++){
        new_size = snprintf(NULL, 0,"%p ", fifo->buf[i]);
        size = new_size > size ? new_size:size;
    }
    string = calloc(fifo->size+1, size);
    for(size_t i=0, size = 0; i<fifo->size; i++)
        size += sprintf(string + size, "%p ", fifo->buf[i]);
    printf("%s\n", string);
    free(string);
    return FIFO_SUCCESS;
}
