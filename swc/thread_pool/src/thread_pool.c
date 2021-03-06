#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "thread_pool.h"
#include "fifo.h"

static void * _slave_func(void *args);

enum thread_state{
    THREAD_STARTING=0,
    THREAD_READY,
    THREAD_BUSY,
    THREAD_STOPPED
};

enum work_flag{
    WORK_STOP=0,
    WORK_WORK,
    WORK_INDEX,
};


static void * _slave_func(void *args)
{
    slave_t *self = args;
    thread_pool_work_t *work;
    int ret;
    int index;
    while(1){
        pthread_mutex_lock(&self->mutex);
        self->state = THREAD_READY;
        pthread_cond_signal(&self->cond);
        pthread_mutex_unlock(&self->mutex);
        work = NULL;
        ret = fifo_pop_index_cond(self->queue, (void **)&work, &index, FIFO_WAIT, &self->state, THREAD_BUSY);
        if(ret != FIFO_SUCCESS){
            fprintf(stderr, "\x1b[31m%s:%s:%d: fifo_pop: error %d\x1b[0m", __FILE__, __func__, __LINE__, ret);
            continue;
        }
        if(!work){
            fprintf(stderr, "\x1b[31m%s:%s:%d: received NULL work\x1b[0m\n", __FILE__, __func__, __LINE__);
            continue;
        }
        if(work->flag == WORK_STOP){
            free(work);
            break;
        } else if(work->flag == WORK_WORK && work->func){
                (work->func)(work->args);
                free(work);
        } else if(work->flag == WORK_INDEX && work->func_index){
                (work->func_index)(work->args, index);
        } else {
            fprintf(stderr, "\x1b[31m%s:%s:%d: fifo_pop: inconsistency (work->flag= %d)\x1b[0m\n", __FILE__, __func__, __LINE__, work->flag);
            continue;
        }
    }
    self->state = THREAD_STOPPED;
    return NULL;
}

int thread_pool_create(thread_pool_t *thread_pool, unsigned int num_slaves, pthread_attr_t *attr)
{
    int ret;
    if(!thread_pool)return THREAD_POOL_UNALLOCATED;
    thread_pool->queue = malloc(sizeof(fifo_t));
    ret = fifo_init(thread_pool->queue, 10*num_slaves);
    if(ret != FIFO_SUCCESS){
        fprintf(stderr, "%s:%s:%d: error %d\n", __FILE__, __func__, __LINE__, ret);
        goto err_queue_open;
    }
    thread_pool->num_slaves = 0;
    thread_pool->slaves = malloc(num_slaves * sizeof(slave_t));
    if(thread_pool->slaves == NULL){
        fprintf(stderr, "%s:%s:%d: ", __FILE__, __func__, __LINE__);
        perror(NULL);
        goto err_slaves;
    }
    for (unsigned int i = 0; i < num_slaves; i++){
        thread_pool->slaves[i].state = THREAD_STARTING;
        thread_pool->slaves[i].queue = thread_pool->queue;
        pthread_mutex_init(&thread_pool->slaves[i].mutex, NULL);
        pthread_cond_init(&thread_pool->slaves[i].cond, NULL);
        ret = pthread_create(&thread_pool->slaves[i].id, attr, _slave_func, &thread_pool->slaves[i]);
        if(ret != 0){
            fprintf(stderr, "%s:%s:%d: error %d\n", __FILE__, __func__, __LINE__, ret);
            goto err_thread;
        }
        thread_pool->num_slaves++;
    }
    goto success;
err_thread:
    for (unsigned int i = 0; i<thread_pool->num_slaves; i++){
        ret = pthread_cancel(thread_pool->slaves[i].id);
        if(ret != 0)
            fprintf(stderr, "\x1b[31m%s:%s:%d: %d\x1b[0m\n", __FILE__, __func__, __LINE__, ret);
    }
    for (unsigned int i = 0; i<thread_pool->num_slaves; i++){
        ret = pthread_join(thread_pool->slaves[i].id, NULL);
        if (ret != 0)
            fprintf(stderr, "\x1b[31m%s:%s:%d: %d\x1b[0m\n", __FILE__, __func__, __LINE__, ret);
    }
    free(thread_pool->slaves);
err_slaves:
    ret = fifo_destroy(thread_pool->queue);
    if(ret != FIFO_SUCCESS)
        fprintf(stderr, "\x1b[31m%s:%s:%d: %d\x1b[0m\n", __FILE__, __func__, __LINE__, ret);
    free(thread_pool->queue);
err_queue_open:
    return THREAD_POOL_KO;
success:
    return THREAD_POOL_OK;
}

int thread_pool_queue(thread_pool_t *thread_pool, void (*func)(void *), void *args)
{
    int ret;
    if(!thread_pool)return THREAD_POOL_UNALLOCATED;
    if(!func)return THREAD_POOL_NULLPTR;
    thread_pool_work_t *work = malloc(sizeof(*work));
    work->flag = WORK_WORK;
    work->func = func;
    work->args = args;
    ret = fifo_push(thread_pool->queue, (void *)work, FIFO_WAIT);
    if(ret != FIFO_SUCCESS){
        fprintf(stderr, "\x1b[31m%s:%s:%d: %d\x1b[0m\n", __FILE__, __func__, __LINE__, ret);
        return THREAD_POOL_KO;
    }
    return THREAD_POOL_OK;
}

int thread_pool_queue_work(thread_pool_t *thread_pool, thread_pool_work_t *work, int index)
{
    int ret;
    if(!thread_pool)return THREAD_POOL_UNALLOCATED;
    work->flag = WORK_INDEX;
    ret = fifo_push_index(thread_pool->queue, (void *)work, index, FIFO_WAIT);
    if(ret != FIFO_SUCCESS){
        fprintf(stderr, "\x1b[31m%s:%s:%d: %d\x1b[0m\n", __FILE__, __func__, __LINE__, ret);
        return THREAD_POOL_KO;
    }
    return THREAD_POOL_OK;
}

int thread_pool_wait(thread_pool_t *thread_pool)
{
    int ret = fifo_wait_empty(thread_pool->queue);
    if(ret != FIFO_SUCCESS){
        fprintf(stderr, "\x1b[31m%s:%s:%d: %d\x1b[0m\n", __FILE__, __func__, __LINE__, ret);
        return THREAD_POOL_KO;
    }
    for (unsigned int i = 0; i<thread_pool->num_slaves; i++){
        pthread_mutex_lock(&thread_pool->slaves[i].mutex);
        if (thread_pool->slaves[i].state != THREAD_READY)
            pthread_cond_wait(&thread_pool->slaves[i].cond, &thread_pool->slaves[i].mutex);
        pthread_mutex_unlock(&thread_pool->slaves[i].mutex);
    }
    return THREAD_POOL_OK;
}

int thread_pool_destroy(thread_pool_t *thread_pool)
{
    int ret;
    thread_pool_work_t *work;
    int current_size = 0;
    for (unsigned int i = 0; i<thread_pool->num_slaves; i++){
        work = malloc(sizeof(thread_pool_work_t));
        work->flag = WORK_STOP;
        ret = fifo_push(thread_pool->queue, (void *)work, FIFO_WAIT);
        if(ret != FIFO_SUCCESS){
            fprintf(stderr, "\x1b[31m%s:%s:%d: %d\x1b[0m\n", __FILE__, __func__, __LINE__, ret);
            return THREAD_POOL_KO;
        }
    }
    for (unsigned int i = 0; i<thread_pool->num_slaves; i++){
       ret = pthread_join(thread_pool->slaves[i].id, NULL); 
       if(ret != 0){
            fprintf(stderr, "\x1b[31m%s:%s:%d: %d\x1b[0m\n", __FILE__, __func__, __LINE__, ret);
            perror(NULL);
            return THREAD_POOL_KO;
        }
        pthread_cond_destroy(&thread_pool->slaves[i].cond);
        pthread_mutex_destroy(&thread_pool->slaves[i].mutex);
    }
    ret = fifo_current_size(thread_pool->queue, &current_size);
    if(ret != FIFO_SUCCESS){
        fprintf(stderr, "\x1b[31m%s:%s:%d: %d\x1b[0m\n", __FILE__, __func__, __LINE__, ret);
        return THREAD_POOL_KO;
    }
    if(current_size > 0)
        fprintf(stderr, "WARNING: Some work will be deleted\n");
    while (fifo_pop(thread_pool->queue, (void **)&work, FIFO_NO_WAIT) != FIFO_EMPTY)
        free(work);
    if(fifo_destroy(thread_pool->queue) != FIFO_SUCCESS){
        fprintf(stderr, "%s:%s:%d: ", __FILE__, __func__, __LINE__);
        return THREAD_POOL_KO;
    }
    free(thread_pool->slaves);
    free(thread_pool->queue);
    return THREAD_POOL_OK;
}
