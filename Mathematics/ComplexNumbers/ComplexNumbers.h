#pragma once

#ifndef COMPLEX_NUMBERS_H
#define COMPLEX_NUMBERS_H

// Complex number structure
typedef struct
{
    float real;
    float imag;
} Complex;

/**
 * Adds two complex numbers.
 * @param a First complex number
 * @param b Second complex number
 * @return The sum of the two complex numbers
 */
Complex complexAdd(Complex *a, Complex *b);

/**
 * Subtracts two complex numbers.
 * @param a First complex number
 * @param b Second complex number
 * @return The difference of the two complex numbers
 */
Complex complexSubtract(Complex *a, Complex *b);

/**
 * Multiplies two complex numbers.
 * @param a First complex number
 * @param b Second complex number
 * @return The product of the two complex numbers
 */
Complex complexMultiply(Complex *a, Complex *b);

/**
 * Divides one complex number by another.
 * @param a The complex number to be divided.
 * @param b The complex number to divide by.
 * @return The quotient of the two complex numbers, or a complex number with NaN values if the denominator is zero.
 */
Complex complexDivide(Complex *a, Complex *b);

/**
 * Returns the complex conjugate of a complex number.
 * @param a The complex number = x + yi
 * @return The conjugate of the complex number = x - yi
 */
Complex complexConjugate(Complex *a);

/**
 * Returns the magnitude of a complex number.
 * @param a The complex number = x + yi
 * @return The magnitude of the complex number = sqrt(x^2 + y^2)
 */
float complexMagnitude(Complex *a);

/**
 * Returns the argument of a complex number, which is the angle in radians between
 * the positive x-axis and the line connecting the origin to the point (x, y) in the complex plane.
 * The function uses atan2 to handle all four quadrants and cases where x = 0.
 *
 * @param a The complex number = x + yi
 * @return The argument of the complex number in radians = atan2(y, x) from -PI to PI
 */
float complexArgument(Complex *a);

/**
 * Prints a complex number to the console in the form x + yi
 * @param c The complex number
 */
void printComplex(Complex *c);

#endif // COMPLEX_NUMBERS_H