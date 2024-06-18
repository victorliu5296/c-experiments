#include "ComplexNumbers.h"
#include <stdio.h>
#include <math.h>  // for NAN
#include <errno.h> // for errno

Complex complexAdd(Complex *a, Complex *b)
{
    Complex result;
    result.real = a->real + b->real;
    result.imag = a->imag + b->imag;
    return result;
}

Complex complexSubtract(Complex *a, Complex *b)
{
    Complex result;
    result.real = a->real - b->real;
    result.imag = a->imag - b->imag;
    return result;
}

Complex complexMultiply(Complex *a, Complex *b)
{
    Complex result;
    result.real = a->real * b->real - a->imag * b->imag;
    result.imag = a->real * b->imag + a->imag * b->real;
    return result;
}

Complex complexDivide(Complex *a, Complex *b) {
    Complex result;
    float denominator = b->real * b->real + b->imag * b->imag;

    if (denominator == 0) {
        printf("Error: Division by zero\n");
        errno = EDOM;  // Set the errno to domain error
        result.real = NAN;
        result.imag = NAN;
    }
    else {
        result.real = (a->real * b->real + a->imag * b->imag) / denominator;
        result.imag = (a->imag * b->real - a->real * b->imag) / denominator;
    }

    return result;
}

Complex complexConjugate(Complex *a)
{
    Complex result;
    result.real = a->real;
    result.imag = -a->imag;
    return result;
}

float complexMagnitude(Complex *a)
{
    return sqrt(a->real * a->real + a->imag * a->imag);
}

float complexArgument(Complex *a)
{
    return atan2(a->imag, a->real);
}

void printComplex(Complex *c)
{
    if (c->imag >= 0)
    {
        printf("%.2f + %.2fi\n", c->real, c->imag);
    }
    else
    {
        printf("%.2f - %.2fi\n", c->real, -c->imag);
    }
}