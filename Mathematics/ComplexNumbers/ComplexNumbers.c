#include "ComplexNumbers.h"
#include <stdio.h>
#include <math.h>  // for NAN
#include <errno.h> // for errno

Complex complex_add(Complex *a, Complex *b)
{
    Complex result;
    result.real = a->real + b->real;
    result.imag = a->imag + b->imag;
    return result;
}

Complex complex_subtract(Complex *a, Complex *b)
{
    Complex result;
    result.real = a->real - b->real;
    result.imag = a->imag - b->imag;
    return result;
}

Complex complex_multiply(Complex *a, Complex *b)
{
    Complex result;
    result.real = a->real * b->real - a->imag * b->imag;
    result.imag = a->real * b->imag + a->imag * b->real;
    return result;
}

Complex complex_divide(Complex *a, Complex *b) {
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

Complex complex_conjugate(Complex *a)
{
    Complex result;
    result.real = a->real;
    result.imag = -a->imag;
    return result;
}

float complex_magnitude(Complex *a)
{
    return sqrt(a->real * a->real + a->imag * a->imag);
}

float complex_argument(Complex *a)
{
    return atan2(a->imag, a->real);
}

void complex_print(Complex *c)
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