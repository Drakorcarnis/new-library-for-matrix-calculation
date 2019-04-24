#ifndef CHECK
#define CHECK
int sanity_check(const void *pointer, const char *function_name);
int square_check(const matrix_t *matrix, const char *function_name);
int symetry_check(const matrix_t *matrix, const char *function_name);
#endif