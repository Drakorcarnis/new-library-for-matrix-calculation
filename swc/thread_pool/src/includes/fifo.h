#ifndef FIFO
#define FIFO
#include <signal.h>
#include <unistd.h>
enum { 
    FIFO_FAIL_UNALLOCATED,
    FIFO_FAIL_MEM,
    FIFO_FAIL_MUTEX,
    FIFO_FAIL_COND,
    FIFO_FULL,
    FIFO_EMPTY,
    FIFO_SUCCESS
};
enum {
    FIFO_WAIT,
    FIFO_NO_WAIT
};
typedef struct
{
    size_t size;
    void ** buf;
    int * index_buf;
    unsigned int push_index;
    unsigned int pop_index;
    unsigned int curr_nb_elt;
    pthread_mutex_t curr_nb_mutex;
    pthread_mutex_t push_mutex;
    pthread_mutex_t pop_mutex;
    pthread_mutex_t push_cond_mutex;
    pthread_mutex_t pop_cond_mutex;
    pthread_mutex_t empty_cond_mutex;
    pthread_cond_t push_cond;
    pthread_cond_t pop_cond;
    pthread_cond_t empty_cond;
} fifo_t;

int fifo_init(fifo_t *fifo, size_t size);
int fifo_push(fifo_t *fifo, void *elt, int wait);
int fifo_push_index(fifo_t *fifo, void *elt, int index, int wait);
int fifo_pop(fifo_t *fifo, void **elt, int wait);
int fifo_pop_index_cond(fifo_t *fifo, void **elt, int *index, int wait, int *cond, int value);
int fifo_wait_empty(fifo_t *fifo);
int fifo_current_size(fifo_t *fifo, int* curr_nb_elt);
int fifo_destroy(fifo_t *fifo);
int fifo_show(fifo_t *fifo);
#endif
