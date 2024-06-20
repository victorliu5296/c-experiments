#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

typedef struct
{
    size_t row_count;
    size_t col_count;
    double *data;
} Matrix;

typedef struct
{
    Matrix inverse;
    bool success;
    const char *error_message;
} InverseResult;

typedef struct
{
    Matrix inverse;
    bool success;
    const char *error_message;
    int iteration_count;
    double final_error;
} IterativeInverseResult;

/**
 * @brief Allocate memory for a matrix of given dimensions.
 * @param rows Number of rows in the matrix.
 * @param cols Number of columns in the matrix.
 * @return A newly allocated Matrix structure.
 */
Matrix matrix_allocate(size_t rows, size_t cols);

/**
 * @brief Create a matrix from a 2D array of values.
 * @param rows Number of rows in the matrix.
 * @param cols Number of columns in the matrix.
 * @param values Pointer to the array of values.
 * @return A newly created Matrix structure.
 */
Matrix matrix_create_from_array(size_t rows, size_t cols, const double *values);

/**
 * @brief Create an identity matrix of given size.
 * @param size The size of the square identity matrix.
 * @return A newly created identity Matrix structure.
 */
Matrix matrix_identity(size_t size);

/**
 * @brief Free the memory allocated for a matrix.
 * @param matrix Pointer to the Matrix structure to be deleted.
 */
void matrix_delete(Matrix *matrix);

/**
 * @brief Set the value of a specific entry in the matrix.
 * @param matrix Pointer to the Matrix structure.
 * @param row Row index of the entry.
 * @param col Column index of the entry.
 * @param value Value to be set.
 * @return 1 if successful, 0 if indices are out of bounds.
 */
int matrix_set_entry(Matrix *matrix, size_t row, size_t col, double value);

/**
 * @brief Get the value of a specific entry in the matrix.
 * @param matrix Pointer to the Matrix structure.
 * @param row Row index of the entry.
 * @param col Column index of the entry.
 * @return The value at the specified position, or NaN if indices are out of bounds.
 */
double matrix_get_entry(const Matrix *matrix, size_t row, size_t col);

/**
 * @brief Print the matrix to the console.
 * @param matrix Pointer to the Matrix structure to be printed.
 */
void matrix_print(const Matrix *matrix);

/**
 * @brief Add two matrices.
 * @param matrix1 Pointer to the first Matrix structure.
 * @param matrix2 Pointer to the second Matrix structure.
 * @return A new Matrix structure containing the sum of the two matrices.
 */
Matrix matrix_add(const Matrix *matrix1, const Matrix *matrix2);

/**
 * @brief Subtract two matrices.
 * @param matrix1 Pointer to the first Matrix structure.
 * @param matrix2 Pointer to the second Matrix structure.
 * @return A new Matrix structure containing the difference of the two matrices.
 */
Matrix matrix_subtract(const Matrix *matrix1, const Matrix *matrix2);

/**
 * @brief Multiply two matrices.
 * @param matrix1 Pointer to the first Matrix structure.
 * @param matrix2 Pointer to the second Matrix structure.
 * @return A new Matrix structure containing the product of the two matrices.
 */
Matrix matrix_multiply(const Matrix *matrix1, const Matrix *matrix2);

/**
 * @brief Transpose a matrix.
 * @param matrix Pointer to the Matrix structure to be transposed.
 * @return A new Matrix structure containing the transposed matrix.
 */
Matrix matrix_transpose(const Matrix *matrix);

/**
 * @brief Scale a matrix by a scalar value.
 * @param matrix Pointer to the Matrix structure to be scaled.
 * @param scalar The scalar value to multiply the matrix by.
 * @return A new Matrix structure containing the scaled matrix.
 */
Matrix matrix_scale(const Matrix *matrix, double scalar);

/**
 * @brief Scale a matrix in place by a scalar value.
 * @param matrix Pointer to the Matrix structure to be scaled.
 * @param scalar The scalar value to multiply the matrix by.
 */
void matrix_scale_inplace(Matrix *matrix, double scalar);

/**
 * @brief Calculate the determinant of a square matrix using LU decomposition.
 * @param matrix Pointer to the Matrix structure.
 * @return The determinant of the matrix, or NAN if the matrix is not square.
 */
double matrix_determinant(const Matrix *matrix);

/**
 * @brief Calculate the trace of a square matrix.
 * @param matrix Pointer to the Matrix structure.
 * @return The trace of the matrix.
 */
double matrix_trace(const Matrix *matrix);

/**
 * @brief Calculate the Frobenius norm of a matrix.
 * @param matrix Pointer to the Matrix structure.
 * @return The Frobenius norm of the matrix.
 */
double matrix_frobenius_norm(const Matrix *matrix);

/**
 * @brief Calculate the 1-norm of a matrix.
 * @param matrix Pointer to the Matrix structure.
 * @return The 1-norm of the matrix.
 */
double matrix_norm_1(const Matrix *matrix);

/**
 * @brief Calculate the infinity norm of a matrix.
 * @param matrix Pointer to the Matrix structure.
 * @return The infinity norm of the matrix.
 */
double matrix_norm_inf(const Matrix *matrix);

/**
 * @brief Calculate the minor of a matrix by removing a specified row and column.
 * @param matrix Pointer to the Matrix structure.
 * @param row Row to remove.
 * @param col Column to remove.
 * @return A new Matrix structure containing the minor.
 */
Matrix matrix_minor(const Matrix *matrix, size_t row, size_t col);

/**
 * @brief Calculate the rank of a matrix.
 * @param matrix Pointer to the Matrix structure.
 * @return The rank of the matrix.
 */
size_t matrix_rank(const Matrix *matrix);

/**
 * @brief Check if a matrix is square.
 * @param matrix Pointer to the Matrix structure.
 * @return 1 if the matrix is square, 0 otherwise.
 */
int matrix_is_square(const Matrix *matrix);

/**
 * @brief Check if a matrix is symmetric.
 * @param matrix Pointer to the Matrix structure.
 * @return 1 if the matrix is symmetric, 0 otherwise.
 */
int matrix_is_symmetric(const Matrix *matrix);

/**
 * @brief Raise a square matrix to a power.
 * @param matrix Pointer to the Matrix structure.
 * @param power The power to raise the matrix to.
 * @return A new Matrix structure containing the result of the matrix raised to the given power.
 */
Matrix matrix_power(const Matrix *matrix, int power);

/**
 * @brief Calculate the inverse of a square matrix using the adjugate method.
 * @param matrix Pointer to the Matrix structure.
 * @return An InverseResult structure containing the inverse matrix and status information.
 */
InverseResult matrix_inverse_adjugate(const Matrix *matrix);

/**
 * @brief Calculate the inverse of a square matrix using Gauss-Jordan elimination.
 * @param matrix Pointer to the Matrix structure.
 * @return An InverseResult structure containing the inverse matrix and status information.
 */
InverseResult matrix_inverse_gauss_jordan(const Matrix *matrix);

/**
 * @brief Calculate the inverse of a square matrix using the Newton-Schulz iteration method.
 * @param matrix Pointer to the Matrix structure.
 * @param max_iterations Maximum number of iterations.
 * @param tolerance Convergence tolerance.
 * @return An IterativeInverseResult structure containing the inverse matrix, status information, and iteration details.
 */
IterativeInverseResult matrix_inverse_newton(const Matrix *matrix, int max_iterations, double tolerance);

/**
 * @brief Swap two rows in a matrix.
 * @param matrix Pointer to the Matrix structure.
 * @param row1 Index of the first row to be swapped.
 * @param row2 Index of the second row to be swapped.
 */
void matrix_swap_rows(Matrix *matrix, size_t row1, size_t row2);

/**
 * @brief Concatenate two matrices horizontally.
 * @param matrix1 Pointer to the first Matrix structure.
 * @param matrix2 Pointer to the second Matrix structure.
 * @return A new Matrix structure containing the horizontally concatenated matrices.
 */
Matrix matrix_concatenate_horizontal(const Matrix *matrix1, const Matrix *matrix2);

/**
 * @brief Concatenate two matrices vertically.
 * @param matrix1 Pointer to the first Matrix structure.
 * @param matrix2 Pointer to the second Matrix structure.
 * @return A new Matrix structure containing the vertically concatenated matrices.
 */
Matrix matrix_concatenate_vertical(const Matrix *matrix1, const Matrix *matrix2);

/**
 * @brief Calculate the maximum absolute difference between two matrices.
 * @param a Pointer to the first Matrix structure.
 * @param b Pointer to the second Matrix structure.
 * @return The maximum absolute difference between the matrices.
 */
double matrix_max_abs_diff(const Matrix *a, const Matrix *b);

#endif // MATRIX_H