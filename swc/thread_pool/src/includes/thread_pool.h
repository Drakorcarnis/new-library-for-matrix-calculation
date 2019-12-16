#ifndef THREAD_POOL
#define THREAD_POOL
#include <pthread.h>
#include "fifo.h"
typedef struct {
    pthread_t       id;
    int             state;
    fifo_t          *queue;
    pthread_mutex_t mutex;
    pthread_cond_t  cond;
} slave_t;

typedef struct {
    slave_t        *slaves;
    unsigned int    num_slaves;
    fifo_t          *queue;
} thread_pool_t;

typedef struct {
    int flag;
    void (*func)(void *);
    void (*func_index)(void *, int);
    void *args;
} thread_pool_work_t;

enum thread_pool_ret{
    THREAD_POOL_KO,
    THREAD_POOL_UNALLOCATED,
    THREAD_POOL_NULLPTR,
    THREAD_POOL_OK
};

int thread_pool_create(thread_pool_t *thread_pool, unsigned int num_slaves, pthread_attr_t *attr);
int thread_pool_queue(thread_pool_t *thread_pool, void (*func)(void *), void *args);
int thread_pool_queue_work(thread_pool_t *thread_pool, thread_pool_work_t *work, int index);
int thread_pool_wait(thread_pool_t *thread_pool);
int thread_pool_destroy(thread_pool_t *thread_pool);
#endif