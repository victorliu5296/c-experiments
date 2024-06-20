#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

Matrix matrix_allocate(size_t rows, size_t cols)
{
    Matrix matrix = {rows, cols, NULL};
    matrix.data = (double *)calloc(rows * cols, sizeof(double));
    if (!matrix.data)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    return matrix;
}

Matrix matrix_create_from_array(size_t rows, size_t cols, const double *values)
{
    Matrix matrix = matrix_allocate(rows, cols);
    for (size_t i = 0; i < rows * cols; ++i)
    {
        matrix.data[i] = values[i];
    }
    return matrix;
}

Matrix matrix_identity(size_t size)
{
    Matrix identity = matrix_allocate(size, size);
    for (size_t i = 0; i < size; ++i)
    {
        identity.data[i * size + i] = 1.0;
    }
    return identity;
}

void matrix_delete(Matrix *matrix)
{
    if (matrix)
    {
        free(matrix->data);
        matrix->data = NULL;
        matrix->row_count = matrix->col_count = 0;
    }
}

int matrix_set_entry(Matrix *matrix, size_t row, size_t col, double value)
{
    if (row >= matrix->row_count || col >= matrix->col_count)
    {
        return 0; // Error: out of bounds
    }
    matrix->data[row * matrix->col_count + col] = value;
    return 1; // Success
}

double matrix_get_entry(const Matrix *matrix, size_t row, size_t col)
{
    if (row >= matrix->row_count || col >= matrix->col_count)
    {
        return NAN; // Error: out of bounds
    }
    return matrix->data[row * matrix->col_count + col];
}

void matrix_print(const Matrix *matrix)
{
    printf("[\n");
    for (size_t i = 0; i < matrix->row_count; ++i)
    {
        printf("  ");
        for (size_t j = 0; j < matrix->col_count; ++j)
        {
            printf("%8.3f", matrix->data[i * matrix->col_count + j]);
            if (j < matrix->col_count - 1)
            {
                printf(", ");
            }
        }
        printf("\n");
    }
    printf("]\n");
}

Matrix matrix_add(const Matrix *matrix1, const Matrix *matrix2)
{
    if (matrix1->row_count != matrix2->row_count || matrix1->col_count != matrix2->col_count)
    {
        fprintf(stderr, "Matrix dimensions do not match for addition.\n");
        return matrix_allocate(0, 0);
    }
    Matrix result = matrix_allocate(matrix1->row_count, matrix1->col_count);
    for (size_t i = 0; i < matrix1->row_count * matrix1->col_count; ++i)
    {
        result.data[i] = matrix1->data[i] + matrix2->data[i];
    }
    return result;
}

Matrix matrix_subtract(const Matrix *matrix1, const Matrix *matrix2)
{
    if (matrix1->row_count != matrix2->row_count || matrix1->col_count != matrix2->col_count)
    {
        fprintf(stderr, "Matrix dimensions do not match for subtraction.\n");
        return matrix_allocate(0, 0);
    }
    Matrix result = matrix_allocate(matrix1->row_count, matrix1->col_count);
    for (size_t i = 0; i < matrix1->row_count * matrix1->col_count; ++i)
    {
        result.data[i] = matrix1->data[i] - matrix2->data[i];
    }
    return result;
}

Matrix matrix_multiply(const Matrix *matrix1, const Matrix *matrix2)
{
    if (matrix1->col_count != matrix2->row_count)
    {
        fprintf(stderr, "Matrix dimensions do not match for multiplication.\n");
        return matrix_allocate(0, 0);
    }
    Matrix result = matrix_allocate(matrix1->row_count, matrix2->col_count);
    for (size_t i = 0; i < matrix1->row_count; ++i)
    {
        for (size_t j = 0; j < matrix2->col_count; ++j)
        {
            double sum = 0.0;
            for (size_t k = 0; k < matrix1->col_count; ++k)
            {
                sum += matrix1->data[i * matrix1->col_count + k] * matrix2->data[k * matrix2->col_count + j];
            }
            result.data[i * result.col_count + j] = sum;
        }
    }
    return result;
}

Matrix matrix_transpose(const Matrix *matrix)
{
    Matrix result = matrix_allocate(matrix->col_count, matrix->row_count);
    for (size_t i = 0; i < matrix->row_count; ++i)
    {
        for (size_t j = 0; j < matrix->col_count; ++j)
        {
            result.data[j * result.col_count + i] = matrix->data[i * matrix->col_count + j];
        }
    }
    return result;
}

Matrix matrix_scale(const Matrix *matrix, double scalar)
{
    Matrix result = matrix_allocate(matrix->row_count, matrix->col_count);
    for (size_t i = 0; i < matrix->row_count * matrix->col_count; ++i)
    {
        result.data[i] = matrix->data[i] * scalar;
    }
    return result;
}

void matrix_scale_inplace(Matrix *matrix, double scalar)
{
    for (size_t i = 0; i < matrix->row_count * matrix->col_count; ++i)
    {
        matrix->data[i] *= scalar;
    }
}

double matrix_determinant(const Matrix *matrix)
{
    if (!matrix_is_square(matrix))
    {
        fprintf(stderr, "Matrix is not square, determinant is not defined.\n");
        return NAN;
    }

    size_t n = matrix->row_count;
    Matrix LU = matrix_allocate(n, n);
    int *P = (int *)malloc(n * sizeof(int));
    if (!P)
    {
        fprintf(stderr, "Memory allocation failed\n");
        matrix_delete(&LU);
        return NAN;
    }

    // Perform LU decomposition with partial pivoting
    for (size_t i = 0; i < n; i++)
    {
        P[i] = i;
        for (size_t j = 0; j < n; j++)
        {
            LU.data[i * n + j] = matrix_get_entry(matrix, i, j);
        }
    }

    int sign = 1;
    for (size_t k = 0; k < n - 1; k++)
    {
        // Find pivot
        size_t pivot_row = k;
        double max_val = fabs(LU.data[k * n + k]);
        for (size_t i = k + 1; i < n; i++)
        {
            double val = fabs(LU.data[i * n + k]);
            if (val > max_val)
            {
                max_val = val;
                pivot_row = i;
            }
        }

        // Swap rows if necessary
        if (pivot_row != k)
        {
            for (size_t j = 0; j < n; j++)
            {
                double temp = LU.data[k * n + j];
                LU.data[k * n + j] = LU.data[pivot_row * n + j];
                LU.data[pivot_row * n + j] = temp;
            }
            int temp = P[k];
            P[k] = P[pivot_row];
            P[pivot_row] = temp;
            sign = -sign;
        }

        // Check for singularity
        if (fabs(LU.data[k * n + k]) < 1e-10)
        {
            matrix_delete(&LU);
            free(P);
            return 0.0; // Matrix is singular
        }

        // Perform elimination
        for (size_t i = k + 1; i < n; i++)
        {
            LU.data[i * n + k] /= LU.data[k * n + k];
            for (size_t j = k + 1; j < n; j++)
            {
                LU.data[i * n + j] -= LU.data[i * n + k] * LU.data[k * n + j];
            }
        }
    }

    // Calculate determinant
    double det = sign;
    for (size_t i = 0; i < n; i++)
    {
        det *= LU.data[i * n + i];
    }

    matrix_delete(&LU);
    free(P);
    return det;
}

double matrix_trace(const Matrix *matrix)
{
    if (matrix->row_count != matrix->col_count)
    {
        fprintf(stderr, "Matrix is not square, trace is not defined.\n");
        return NAN;
    }
    double trace = 0.0;
    for (size_t i = 0; i < matrix->row_count; ++i)
    {
        trace += matrix->data[i * matrix->col_count + i];
    }
    return trace;
}

double matrix_frobenius_norm(const Matrix *matrix)
{
    double norm = 0.0;
    for (size_t i = 0; i < matrix->row_count * matrix->col_count; ++i)
    {
        norm += matrix->data[i] * matrix->data[i];
    }
    return sqrt(norm);
}

double matrix_norm_1(const Matrix *matrix)
{
    double max_sum = 0.0;
    for (size_t j = 0; j < matrix->col_count; j++)
    {
        double sum = 0.0;
        for (size_t i = 0; i < matrix->row_count; i++)
        {
            sum += fabs(matrix_get_entry(matrix, i, j));
        }
        if (sum > max_sum)
        {
            max_sum = sum;
        }
    }
    return max_sum;
}

double matrix_norm_inf(const Matrix *matrix)
{
    double max_sum = 0.0;
    for (size_t i = 0; i < matrix->row_count; i++)
    {
        double sum = 0.0;
        for (size_t j = 0; j < matrix->col_count; j++)
        {
            sum += fabs(matrix_get_entry(matrix, i, j));
        }
        if (sum > max_sum)
        {
            max_sum = sum;
        }
    }
    return max_sum;
}

Matrix matrix_minor(const Matrix *matrix, size_t row, size_t col)
{
    if (!matrix_is_square(matrix) || row >= matrix->row_count || col >= matrix->col_count)
    {
        // Return an empty matrix if input is invalid
        return matrix_allocate(0, 0);
    }
    size_t n = matrix->row_count;
    Matrix minor = matrix_allocate(n - 1, n - 1);
    size_t minor_i = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (i == row)
            continue;
        size_t minor_j = 0;
        for (size_t j = 0; j < n; ++j)
        {
            if (j == col)
                continue;
            double value = matrix_get_entry(matrix, i, j);
            matrix_set_entry(&minor, minor_i, minor_j, value);
            minor_j++;
        }
        minor_i++;
    }
    return minor;
}

size_t matrix_rank(const Matrix *matrix)
{
    Matrix temp = matrix_allocate(matrix->row_count, matrix->col_count);
    for (size_t i = 0; i < matrix->row_count * matrix->col_count; ++i)
    {
        temp.data[i] = matrix->data[i];
    }

    size_t rank = matrix->col_count;

    for (size_t row = 0; row < rank; ++row)
    {
        if (fabs(temp.data[row * temp.col_count + row]) > 1e-10)
        {
            for (size_t col = 0; col < temp.row_count; ++col)
            {
                if (col != row)
                {
                    double mult = temp.data[col * temp.col_count + row] / temp.data[row * temp.col_count + row];
                    for (size_t i = 0; i < rank; ++i)
                    {
                        temp.data[col * temp.col_count + i] -= mult * temp.data[row * temp.col_count + i];
                    }
                }
            }
        }
        else
        {
            // If the diagonal element is zero, try to swap with a lower row
            int swap_row = -1;
            for (size_t i = row + 1; i < temp.row_count; ++i)
            {
                if (fabs(temp.data[i * temp.col_count + row]) > 1e-10)
                {
                    swap_row = i;
                    break;
                }
            }

            if (swap_row != -1)
            {
                matrix_swap_rows(&temp, row, swap_row);
            }
            else
            {
                // If we couldn't find a non-zero element, reduce the rank
                rank--;
                for (size_t i = 0; i < temp.row_count; ++i)
                {
                    temp.data[i * temp.col_count + row] = temp.data[i * temp.col_count + rank];
                }
                row--;
            }
        }
    }

    matrix_delete(&temp);
    return rank;
}

int matrix_is_square(const Matrix *matrix)
{
    return matrix->row_count == matrix->col_count;
}

int matrix_is_symmetric(const Matrix *matrix)
{
    if (!matrix_is_square(matrix))
    {
        return 0;
    }

    for (size_t i = 0; i < matrix->row_count; ++i)
    {
        for (size_t j = i + 1; j < matrix->col_count; ++j)
        {
            if (fabs(matrix->data[i * matrix->col_count + j] - matrix->data[j * matrix->col_count + i]) > 1e-10)
            {
                return 0;
            }
        }
    }
    return 1;
}

Matrix matrix_power(const Matrix *matrix, int power)
{
    if (!matrix_is_square(matrix))
    {
        fprintf(stderr, "Matrix must be square for power operation.\n");
        return matrix_allocate(0, 0);
    }

    if (power < 0)
    {
        fprintf(stderr, "Negative power is not supported.\n");
        return matrix_allocate(0, 0);
    }

    if (power == 0)
    {
        return matrix_identity(matrix->row_count);
    }

    Matrix result = matrix_allocate(matrix->row_count, matrix->col_count);
    for (size_t i = 0; i < matrix->row_count * matrix->col_count; ++i)
    {
        result.data[i] = matrix->data[i];
    }

    for (int i = 1; i < power; ++i)
    {
        Matrix temp = matrix_multiply(&result, matrix);
        matrix_delete(&result);
        result = temp;
    }

    return result;
}

InverseResult matrix_inverse_adjugate(const Matrix *matrix)
{
    InverseResult result = {matrix_allocate(0, 0), false, NULL};
    if (!matrix_is_square(matrix))
    {
        result.error_message = "Matrix is not square, inverse is not defined.";
        return result;
    }
    double det = matrix_determinant(matrix);
    if (fabs(det) < 1e-10)
    {
        result.error_message = "Matrix is singular, inverse is not defined.";
        return result;
    }
    size_t n = matrix->row_count;
    Matrix adjoint = matrix_allocate(n, n);
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            Matrix submatrix = matrix_minor(matrix, i, j);
            double cofactor = ((i + j) % 2 == 0 ? 1 : -1) * matrix_determinant(&submatrix);
            matrix_set_entry(&adjoint, j, i, cofactor);
            matrix_delete(&submatrix);
        }
    }
    result.inverse = matrix_scale(&adjoint, 1.0 / det);
    matrix_delete(&adjoint);
    result.success = true;
    return result;
}

InverseResult matrix_inverse_gauss_jordan(const Matrix *matrix)
{
    InverseResult result = {matrix_allocate(0, 0), false, NULL};
    if (!matrix_is_square(matrix))
    {
        result.error_message = "Matrix is not square, inverse is not defined.";
        return result;
    }

    size_t n = matrix->row_count;
    Matrix identity = matrix_identity(n);
    Matrix augmented = matrix_concatenate_horizontal(matrix, &identity);

    for (size_t i = 0; i < n; ++i)
    {
        // Find pivot
        size_t pivot_row = i;
        for (size_t j = i + 1; j < n; ++j)
        {
            if (fabs(matrix_get_entry(&augmented, j, i)) > fabs(matrix_get_entry(&augmented, pivot_row, i)))
            {
                pivot_row = j;
            }
        }

        // Swap rows if necessary
        if (pivot_row != i)
        {
            matrix_swap_rows(&augmented, i, pivot_row);
        }

        // Check for singularity
        double pivot = matrix_get_entry(&augmented, i, i);
        if (fabs(pivot) < 1e-10)
        {
            matrix_delete(&augmented);
            result.error_message = "Matrix is singular and cannot be inverted.";
            return result;
        }

        // Normalize pivot row
        for (size_t j = 0; j < 2 * n; ++j)
        {
            double value = matrix_get_entry(&augmented, i, j) / pivot;
            matrix_set_entry(&augmented, i, j, value);
        }

        // Eliminate pivot column
        for (size_t k = 0; k < n; ++k)
        {
            if (k != i)
            {
                double factor = matrix_get_entry(&augmented, k, i);
                for (size_t j = 0; j < 2 * n; ++j)
                {
                    double value = matrix_get_entry(&augmented, k, j) - factor * matrix_get_entry(&augmented, i, j);
                    matrix_set_entry(&augmented, k, j, value);
                }
            }
        }
    }

    // Extract inverse from augmented matrix
    result.inverse = matrix_allocate(n, n);
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            double value = matrix_get_entry(&augmented, i, j + n);
            matrix_set_entry(&result.inverse, i, j, value);
        }
    }

    matrix_delete(&augmented);
    matrix_delete(&identity);
    result.success = true;
    return result;
}

IterativeInverseResult matrix_inverse_newton(const Matrix *matrix, int max_iterations, double tolerance)
{
    IterativeInverseResult result;
    result.success = false;
    result.error_message = "Unknown error";
    result.iteration_count = 0;
    result.final_error = 0.0;

    if (!matrix_is_square(matrix))
    {
        result.error_message = "Matrix is not square";
        return result;
    }

    size_t n = matrix->row_count;
    Matrix X = matrix_transpose(matrix);
    double frobenius_norm_sq = matrix_frobenius_norm(matrix) * matrix_frobenius_norm(matrix);
    matrix_scale_inplace(&X, 1.0 / frobenius_norm_sq);

    Matrix I = matrix_identity(n);

    Matrix A = *matrix;
    Matrix AX = matrix_multiply(&A, &X);
    Matrix I_minus_AX = matrix_subtract(&I, &AX);
    double error = matrix_max_abs_diff(&I, &AX);

    for (int k = 0; k < max_iterations; ++k)
    {
        if (error < tolerance)
        {
            result.success = true;
            result.error_message = "Convergence achieved";
            break;
        }

        Matrix correction = matrix_multiply(&X, &I_minus_AX);
        Matrix X_new = matrix_add(&X, &correction);

        matrix_delete(&correction);
        matrix_delete(&AX);
        matrix_delete(&I_minus_AX);

        AX = matrix_multiply(&A, &X_new);
        I_minus_AX = matrix_subtract(&I, &AX);
        error = matrix_max_abs_diff(&I, &AX);

        matrix_delete(&X);
        X = X_new;

        result.iteration_count = k + 1;
        result.final_error = error;
    }

    matrix_delete(&AX);
    matrix_delete(&I_minus_AX);
    matrix_delete(&I);

    result.inverse = X;
    if (!result.success)
    {
        result.error_message = "Failed to converge within the maximum number of iterations";
    }

    return result;
}

void matrix_swap_rows(Matrix *matrix, size_t row1, size_t row2)
{
    if (row1 >= matrix->row_count || row2 >= matrix->row_count)
    {
        fprintf(stderr, "Invalid row indices for swapping.\n");
        return;
    }

    for (size_t i = 0; i < matrix->col_count; ++i)
    {
        double temp = matrix->data[row1 * matrix->col_count + i];
        matrix->data[row1 * matrix->col_count + i] = matrix->data[row2 * matrix->col_count + i];
        matrix->data[row2 * matrix->col_count + i] = temp;
    }
}

Matrix matrix_concatenate_horizontal(const Matrix *matrix1, const Matrix *matrix2)
{
    if (matrix1->row_count != matrix2->row_count)
    {
        fprintf(stderr, "Matrices must have the same number of rows for horizontal concatenation.\n");
        return matrix_allocate(0, 0);
    }

    Matrix result = matrix_allocate(matrix1->row_count, matrix1->col_count + matrix2->col_count);

    for (size_t i = 0; i < matrix1->row_count; ++i)
    {
        for (size_t j = 0; j < matrix1->col_count; ++j)
        {
            result.data[i * result.col_count + j] = matrix1->data[i * matrix1->col_count + j];
        }
        for (size_t j = 0; j < matrix2->col_count; ++j)
        {
            result.data[i * result.col_count + matrix1->col_count + j] = matrix2->data[i * matrix2->col_count + j];
        }
    }

    return result;
}

Matrix matrix_concatenate_vertical(const Matrix *matrix1, const Matrix *matrix2)
{
    if (matrix1->col_count != matrix2->col_count)
    {
        fprintf(stderr, "Matrices must have the same number of columns for vertical concatenation.\n");
        return matrix_allocate(0, 0);
    }

    Matrix result = matrix_allocate(matrix1->row_count + matrix2->row_count, matrix1->col_count);

    for (size_t i = 0; i < matrix1->row_count; ++i)
    {
        for (size_t j = 0; j < matrix1->col_count; ++j)
        {
            result.data[i * result.col_count + j] = matrix1->data[i * matrix1->col_count + j];
        }
    }

    for (size_t i = 0; i < matrix2->row_count; ++i)
    {
        for (size_t j = 0; j < matrix2->col_count; ++j)
        {
            result.data[(i + matrix1->row_count) * result.col_count + j] = matrix2->data[i * matrix2->col_count + j];
        }
    }

    return result;
}

double matrix_max_abs_diff(const Matrix *matrix1, const Matrix *matrix2)
{
    if (matrix1->row_count != matrix2->row_count || matrix1->col_count != matrix2->col_count)
    {
        fprintf(stderr, "Matrices must have the same dimensions for max_abs_diff.\n");
        return NAN;
    }

    double max_diff = 0.0;
    for (size_t i = 0; i < matrix1->row_count * matrix1->col_count; ++i)
    {
        double diff = fabs(matrix1->data[i] - matrix2->data[i]);
        if (diff > max_diff)
        {
            max_diff = diff;
        }
    }

    return max_diff;
}