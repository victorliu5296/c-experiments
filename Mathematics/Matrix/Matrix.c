#include <stdio.h>
#include <stdlib.h>
#include "Matrix.h"

// Function to create a matrix
Matrix createMatrix(int rows, int cols)
{
    Matrix matrix;
    matrix.rows = rows;
    matrix.cols = cols;
    matrix.data = (double **)malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++)
    {
        matrix.data[i] = (double *)malloc(cols * sizeof(double));
    }
    return matrix;
}

// Function to create a matrix from a 2D array
Matrix createMatrixFromArray(int rows, int cols, double *values)
{
    Matrix matrix = createMatrix(rows, cols);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            matrix.data[i][j] = values[i * cols + j];
        }
    }
    return matrix;
}

// Function to free the memory allocated for a matrix
void freeMatrix(Matrix *matrix)
{
    for (int i = 0; i < matrix->rows; i++)
    {
        free(matrix->data[i]);
    }
    free(matrix->data);
    matrix->data = NULL;
    matrix->rows = 0;
    matrix->cols = 0;
}

// Function to set a value in the matrix
void setMatrixValue(Matrix *matrix, int row, int col, double value)
{
    if (row >= 0 && row < matrix->rows && col >= 0 && col < matrix->cols)
    {
        matrix->data[row][col] = value;
    }
}

// Function to get a value from the matrix
double getMatrixValue(Matrix *matrix, int row, int col)
{
    if (row >= 0 && row < matrix->rows && col >= 0 && col < matrix->cols)
    {
        return matrix->data[row][col];
    }
    return 0.0; // Return 0.0 if the indices are out of bounds
}

// Function to print the matrix
void printMatrix(Matrix *matrix)
{
    for (int i = 0; i < matrix->rows; i++)
    {
        for (int j = 0; j < matrix->cols; j++)
        {
            printf("%lf ", matrix->data[i][j]);
        }
        printf("\n");
    }
}