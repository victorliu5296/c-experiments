#include <stdio.h>
#include "Matrix.h"
#include "examples.h"

void matrix()
{
    // Define a 2D array of values
    double values[3][3] = 
    {
        {1.0, 2.0, 3.0},
        {4.0, 5.0, 6.0},
        {7.0, 8.0, 9.0}
    };

    // Create a matrix from the 2D array
    Matrix matrix = createMatrixFromArray(3, 3, (double*)values);

    // Print the matrix
    printMatrix(&matrix);

    // Free the memory allocated for the matrix
    freeMatrix(&matrix);
}