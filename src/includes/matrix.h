#ifndef MATRIX
#define MATRIX

typedef struct {
    int rows;
    int columns;
    double **coeff;
} matrix_t;

// PLU factorisation 
typedef struct {
    int nb_perm;
    matrix_t *P;
    matrix_t *L;
    matrix_t *U;
} plu_t;

// Cholesky factorisation output
typedef struct {
    int nb_perm;
    matrix_t *L;
    matrix_t *D;
} cholesky_t;

// Matrix creation functions
matrix_t *  matrix_create(int rows, int columns);
matrix_t *  matrix_identity(int n);
matrix_t *  matrix_permutation(int line1, int line2, int n);                                // Return a permutation matrix of rank n for two lines
plu_t *     plu_create(int rank);
matrix_t *  matrix_copy(const matrix_t *matrix);
matrix_t *  str2matrix(int argc, char **argv, char separator);
matrix_t *  file2matrix(char *filename);

// Matrix destruction functions
void        matrix_free(matrix_t *matrix);
void        plu_free(plu_t *plu);

// Matrix display function
void        matrix_display(const matrix_t *matrix);                                         // Display matrix representation to stdout

// Matrix computation functions

// Basic operations
matrix_t *  matrix_transp_f(const matrix_t *matrix);                                        // Return transposed matrix
matrix_t *  matrix_add_f(const matrix_t *matrix1, const matrix_t *matrix2);                 // Return matrix1 + matrix2
matrix_t *  matrix_mult_scalar_f(const matrix_t *matrix, double lambda);                    // Return λ * matrix
matrix_t *  matrix_mult_f(const matrix_t *matrix1, const matrix_t *matrix2);                // Return matrix1 * matrix2
matrix_t *  matrix_pow_f(const matrix_t *matrix, int pow);                                  // Return matrix^pow
matrix_t *  matrix_shrink_f(const matrix_t *matrix, int skipped_row, int skipped_column);   // Return matrix with cut skipped row and skipped column
matrix_t *  matrix_solve_diag_inf(const matrix_t *A, const matrix_t *B);                    // Resolve AX=B with A being an inferior-diagonal matrix. Return X
matrix_t *  matrix_solve_diag_sup(const matrix_t *A, const matrix_t *B);                    // Resolve AX=B with A being an superior-diagonal matrix. Return X

// Raw methods. For fun only. Do never use them, cuz you've NO reason to use them. Really.
double      matrix_det_raw_f(const matrix_t *matrix);                                       // Return |matrix| with brute force method
matrix_t *  matrix_inverse_raw_f(const matrix_t *matrix);                                   // Return matrix^-1 computed with brute force method
matrix_t *  matrix_solve_raw_f(const matrix_t *A, const matrix_t *B);                       // Resolve AX=B with brute force method. Return X
matrix_t *  matrix_com_f(const matrix_t *matrix);                                           // Return comatrix
matrix_t *  matrix_comp_f(const matrix_t *matrix);                                          // Return complementary matrix

// Highly optimized fast methods based upon PLU decomposition for square matrix.
plu_t *     matrix_plu_f(const matrix_t *matrix);                                           // Return PLU decomposition of matrix
double      matrix_det_plu_f(const matrix_t *matrix);                                       // Return |matrix| with PLU method
matrix_t *  matrix_solve_plu_f(const matrix_t *A, const matrix_t *B);                       // Resolve AX=B with PLU method. Return X
matrix_t *  matrix_inverse_plu_f(const matrix_t *matrix);                                   // Return matrix^-1 computed with PLU method

// Highly optimized fast methods based upon Cholesky decomposition for square symetric matrix.
matrix_t * matrix_cholesky_f(const matrix_t *matrix);                                       // Return L with A=L*transp(L) using Cholesky factorisation
matrix_t * matrix_solve_cholesky_f(const matrix_t *A, const matrix_t *B);                   // Resolve AX=B with Cholesky method. Return X
matrix_t * matrix_inverse_cholesky_f(const matrix_t *matrix);                               // Return matrix^-1 computed with Cholesky method
#endif