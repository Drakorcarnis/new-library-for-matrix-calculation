#ifndef MATRIX
#define MATRIX
#include <complex.h>
typedef struct {
    int rows;
    int columns;
    double complex **coeff;
} matrix_t;


// Matrix creation functions
matrix_t *  matrix_create(int rows, int columns);                                           // Creates a 0-filled rows*columns matrix
matrix_t *  matrix_symetric_random(int rows, int columns);                                  // Creates a 0-filled rows*columns symetric matrix
matrix_t *  matrix_random(int rows, int columns);                                           // Creates a random-filled rows*columns matrix (-10<Aij<10)
matrix_t *  matrix_identity(int n);                                                         // Creates Identity matrix of rank n
matrix_t *  matrix_permutation(int line1, int line2, int n);                                // Creates a permutation matrix of rank n for two lines
matrix_t *  matrix_copy(const matrix_t *matrix);                                            // Copies a matrix
matrix_t *  str2matrix(int argc, char **argv, char separator);                              // Creates a matrix from 'separator' separated numbers from 'argc' strings
matrix_t *  file2matrix(char *filename);                                                    // Creates a matrix from a file. Lines separated by line-feed, numbers separated by spaces

// Matrix destruction functions
void        matrix_free(matrix_t *matrix);                                                  // Destroys a matrix

// Matrix display function
void        matrix_display(const matrix_t *matrix);                                         // Display matrix representation to stdout with standard precision
void        matrix_display_exact(const matrix_t *matrix, int precision);                    // Display matrix representation to stdout with specified precision

// Matrix computation functions

// Basic operations
matrix_t *  matrix_transp_f(const matrix_t *matrix);                                        // Return transposed matrix
matrix_t *  matrix_add_f(const matrix_t *matrix1, const matrix_t *matrix2);                 // Return matrix1 + matrix2
matrix_t *  matrix_mult_scalar_f(const matrix_t *matrix, double complex lambda);            // Return λ * matrix
matrix_t *  matrix_mult_f(const matrix_t *matrix1, const matrix_t *matrix2);                // Return matrix1 * matrix2
matrix_t *  matrix_pow_f(const matrix_t *matrix, int pow);                                  // Return matrix^pow
matrix_t *  matrix_shrink_f(const matrix_t *matrix, int skipped_row, int skipped_column);   // Return copy of matrix without 'skipped_row' and 'skipped_column'
matrix_t *  matrix_solve_diag_inf(const matrix_t *A, const matrix_t *B);                    // Resolve AX=B with A being an inferior-diagonal matrix. Return X
matrix_t *  matrix_solve_diag_sup(const matrix_t *A, const matrix_t *B);                    // Resolve AX=B with A being an superior-diagonal matrix. Return X

// Raw methods. For fun only. Do never use them, cuz you've NO reason to use them. Really.
double complex matrix_det_raw_f(const matrix_t *matrix);                               // Return |matrix| with brute force method
matrix_t *  matrix_inverse_raw_f(const matrix_t *matrix);                                   // Return matrix^-1 computed with brute force method
matrix_t *  matrix_solve_raw_f(const matrix_t *A, const matrix_t *B);                       // Resolve AX=B with brute force method. Return X
matrix_t *  matrix_com_f(const matrix_t *matrix);                                           // Return comatrix
matrix_t *  matrix_comp_f(const matrix_t *matrix);                                          // Return complementary matrix

// Highly optimized fast methods based upon PLU decomposition for square matrix.
double complex matrix_det_plu_f(const matrix_t *matrix);                               // Return |matrix| with PLU method
matrix_t *  matrix_solve_plu_f(const matrix_t *A, const matrix_t *B);                       // Resolve AX=B with PLU method. Return X
matrix_t *  matrix_inverse_plu_f(const matrix_t *matrix);                                   // Return matrix^-1 computed with PLU method

// Highly optimized fast methods based upon Cholesky decomposition for square symetric matrix.
double complex matrix_det_cholesky_f(const matrix_t *matrix);
matrix_t * matrix_solve_cholesky_f(const matrix_t *A, const matrix_t *B);                   // Resolve AX = B with Cholesky method. Return X
matrix_t * matrix_inverse_cholesky_f(const matrix_t *matrix);                               // Return matrix^-1 computed with Cholesky method
#endif