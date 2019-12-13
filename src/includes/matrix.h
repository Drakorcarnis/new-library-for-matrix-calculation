#ifndef MATRIX
#define MATRIX
#define TYPE double

typedef struct {
    size_t rows;
    size_t columns;
    TYPE **coeff;
} matrix_t;

// Library initialisation
int         libmatrix_init(void);
int         libmatrix_end(void);
// Matrix creation functions
matrix_t *  matrix_create(size_t rows, size_t columns);                         // Creates a 0-filled rows*columns matrix
matrix_t *  matrix_identity(size_t n);                                                // Creates Identity matrix of rank n
matrix_t *  matrix_permutation(size_t line1, size_t line2, size_t n);     // Creates a permutation matrix of rank n for two lines
matrix_t *  matrix_copy(const matrix_t *matrix);                                            // Copies a matrix

// Matrix destruction functions
void        matrix_free(matrix_t *matrix);                                                  // Destroys a matrix

// Matrix computation functions

// Basic operations
matrix_t *  matrix_transp_f(const matrix_t *matrix);                                        // Return transposed matrix
matrix_t *  matrix_add_f(const matrix_t *matrix1, const matrix_t *matrix2);                 // Return matrix1 + matrix2
matrix_t *  matrix_mult_scalar_f(const matrix_t *matrix, TYPE lambda);                      // Return λ * matrix
matrix_t *  matrix_mult_f(const matrix_t *matrix1, const matrix_t *matrix2);                // Return matrix1 * matrix2
matrix_t *  BAKmatrix_mult_f(const matrix_t *matrix1, const matrix_t *matrix2);                // Return matrix1 * matrix2
matrix_t *  matrix_pow_f(const matrix_t *matrix, int pow);                                  // Return matrix^pow

// Raw methods. For fun only. Do never use them, cuz you've NO reason to use them. Really.
TYPE        matrix_det_raw_f(const matrix_t *matrix);                                       // Return |matrix| with brute force method
matrix_t *  matrix_inverse_raw_f(const matrix_t *matrix);                                   // Return matrix^-1 computed with brute force method
matrix_t *  matrix_solve_raw_f(const matrix_t *A, const matrix_t *B);                       // Resolve AX=B with brute force method. Return X
matrix_t *  matrix_com_f(const matrix_t *matrix);                                           // Return comatrix
matrix_t *  matrix_comp_f(const matrix_t *matrix);                                          // Return complementary matrix

// Highly optimized fast methods based upon PLU decomposition for square matrix.
TYPE        matrix_det_plu_f(const matrix_t *matrix);                                       // Return |matrix| with PLU method
matrix_t *  matrix_solve_plu_f(const matrix_t *A, const matrix_t *B);                       // Resolve AX=B with PLU method. Return X
matrix_t *  matrix_inverse_plu_f(const matrix_t *matrix);                                   // Return matrix^-1 computed with PLU method

// Highly optimized fast methods based upon Cholesky decomposition for square symetric matrix.
TYPE matrix_det_cholesky_f(const matrix_t *matrix);
matrix_t * matrix_solve_cholesky_f(const matrix_t *A, const matrix_t *B);                   // Resolve AX = B with Cholesky method. Return X
matrix_t * matrix_inverse_cholesky_f(const matrix_t *matrix);                               // Return matrix^-1 computed with Cholesky method
#endif