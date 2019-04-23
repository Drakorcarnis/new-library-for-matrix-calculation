#ifndef MATRIX_TOOLS
#define MATRIX_TOOLS
matrix_t *  matrix_random(int rows, int columns);                                           // Creates a random-filled rows*columns matrix (-10<Aij<10)
matrix_t *  matrix_symetric_random(int rows, int columns);                                  // Creates a 0-filled rows*columns symetric matrix
matrix_t *  str2matrix(int argc, char **argv, char separator);                              // Creates a matrix from 'separator' separated numbers from 'argc' strings
matrix_t *  file2matrix(char *filename);                                                    // Creates a matrix from a file. Lines separated by line-feed, numbers separated by spaces

// Matrix display function
void        matrix_display(const matrix_t *matrix);                                         // Display matrix representation to stdout with standard precision
void        matrix_display_exact(const matrix_t *matrix, int precision);                    // Display matrix representation to stdout with specified precision
int         matrix2file(matrix_t *matrix, char * filename);

char * format_time(const long long input_time, char* format);
long long mstime(void);
int test_matrix_equality(const matrix_t *matrix1, const matrix_t *matrix2, int precision);
#endif
