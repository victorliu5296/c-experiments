#include <stdio.h>
#include "matrix.h"

void matrix()
{
    printf("Matrix example running\n\n");

    // Define a 2D array of values
    double values1[3][3] = {
        {4.0, 7.0, 2.0},
        {2.0, 6.0, 3.0},
        {1.0, 8.0, 5.0}};

    double values2[3][3] = {
        {9.0, 8.0, 7.0},
        {6.0, 5.0, 4.0},
        {3.0, 2.0, 1.0}};

    // Create matrices from the 2D arrays
    Matrix matrix1 = matrix_create_from_array(3, 3, (double *)values1);
    Matrix matrix2 = matrix_create_from_array(3, 3, (double *)values2);

    printf("Matrix 1:\n");
    matrix_print(&matrix1);

    printf("\nMatrix 2:\n");
    matrix_print(&matrix2);

    // Addition
    Matrix sum = matrix_add(&matrix1, &matrix2);
    printf("\nSum of matrices:\n");
    matrix_print(&sum);

    // Subtraction
    Matrix difference = matrix_subtract(&matrix1, &matrix2);
    printf("\nDifference of matrices:\n");
    matrix_print(&difference);

    // Multiplication
    Matrix product = matrix_multiply(&matrix1, &matrix2);
    printf("\nProduct of matrices:\n");
    matrix_print(&product);

    // Transpose
    Matrix transpose = matrix_transpose(&matrix1);
    printf("\nTranspose of Matrix 1:\n");
    matrix_print(&transpose);

    // Scaling
    Matrix scaled = matrix_scale(&matrix1, 2.0);
    printf("\nMatrix 1 scaled by 2:\n");
    matrix_print(&scaled);

    // Determinant
    double det = matrix_determinant(&matrix1);
    printf("\nDeterminant of Matrix 1: %f\n", det);

    // Trace
    double trace = matrix_trace(&matrix1);
    printf("\nTrace of Matrix 1: %f\n", trace);

    // Frobenius norm
    double norm = matrix_frobenius_norm(&matrix1);
    printf("\nFrobenius norm of Matrix 1: %f\n", norm);

    // 1-norm
    double norm_1 = matrix_norm_1(&matrix1);
    printf("\n1-norm of Matrix 1: %f\n", norm_1);

    // Infinity norm
    double norm_inf = matrix_norm_inf(&matrix1);
    printf("\nInfinity norm of Matrix 1: %f\n", norm_inf);

    // Inverse (using Adjugate method)
    printf("\nInverse of Matrix 1 (Adjugate method):\n");
    InverseResult inv_result1 = matrix_inverse_adjugate(&matrix1);
    if (inv_result1.success)
    {
        matrix_print(&inv_result1.inverse);
        matrix_delete(&inv_result1.inverse);
    }
    else
    {
        printf("Error: %s\n", inv_result1.error_message);
    }

    // Inverse (using Gauss-Jordan method)
    printf("\nInverse of Matrix 1 (Gauss-Jordan method):\n");
    InverseResult inv_result2 = matrix_inverse_gauss_jordan(&matrix1);
    if (inv_result2.success)
    {
        matrix_print(&inv_result2.inverse);
        matrix_delete(&inv_result2.inverse);
    }
    else
    {
        printf("Error: %s\n", inv_result2.error_message);
    }

    // Inverse (using Newton-Schulz iteration method)
    printf("\nInverse of Matrix 1 (Newton-Schulz method):\n");
    IterativeInverseResult iter_inv_result = matrix_inverse_newton(&matrix1, 100, 1e-6);
    if (iter_inv_result.success)
    {
        matrix_print(&iter_inv_result.inverse);
        printf("Iterations: %d, Final error: %e\n", iter_inv_result.iteration_count, iter_inv_result.final_error);
        matrix_delete(&iter_inv_result.inverse);
    }
    else
    {
        printf("Error: %s\n", iter_inv_result.error_message);
    }

    // Rank
    size_t rank = matrix_rank(&matrix1);
    printf("\nRank of Matrix 1: %zu\n", rank);

    // Identity matrix
    Matrix identity = matrix_identity(3);
    printf("\n3x3 Identity matrix:\n");
    matrix_print(&identity);

    // Check if matrix is square
    printf("\nIs Matrix 1 square? %s\n", matrix_is_square(&matrix1) ? "Yes" : "No");

    // Check if matrix is symmetric
    printf("Is Matrix 1 symmetric? %s\n", matrix_is_symmetric(&matrix1) ? "Yes" : "No");

    // Matrix power
    Matrix power = matrix_power(&matrix1, 2);
    printf("\nMatrix 1 squared:\n");
    matrix_print(&power);

    // Concatenate matrices
    Matrix concat_h = matrix_concatenate_horizontal(&matrix1, &matrix2);
    printf("\nHorizontal concatenation of Matrix 1 and Matrix 2:\n");
    matrix_print(&concat_h);

    Matrix concat_v = matrix_concatenate_vertical(&matrix1, &matrix2);
    printf("\nVertical concatenation of Matrix 1 and Matrix 2:\n");
    matrix_print(&concat_v);

    // Free the memory allocated for all matrices
    matrix_delete(&matrix1);
    matrix_delete(&matrix2);
    matrix_delete(&sum);
    matrix_delete(&difference);
    matrix_delete(&product);
    matrix_delete(&transpose);
    matrix_delete(&scaled);
    matrix_delete(&identity);
    matrix_delete(&power);
    matrix_delete(&concat_h);
    matrix_delete(&concat_v);
}