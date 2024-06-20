#ifndef MATRIX_H
#define MATRIX_H

typedef struct
{
    int rows;
    int cols;
    double **data;
} Matrix;

// Function to create a matrix
Matrix createMatrix(int rows, int cols);

// Function to create a matrix from a 2D array
Matrix createMatrixFromArray(int rows, int cols, double *values);

// Function to free the memory allocated for a matrix
void freeMatrix(Matrix *matrix);

// Function to set a value in the matrix
void setMatrixValue(Matrix *matrix, int row, int col, double value);

// Function to get a value from the matrix
double getMatrixValue(Matrix *matrix, int row, int col);

// Function to print the matrix
void printMatrix(Matrix *matrix);

#endif // MATRIX_H