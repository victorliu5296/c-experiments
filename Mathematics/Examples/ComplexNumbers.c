#include "ComplexNumbers.h"

int main()
{
    Complex a = {.real = 1, .imag = 2};
    Complex b = {.real = 3, .imag = 4};
    Complex c = complexAdd(&a, &b);
    printComplex(&c);

    return 0;
}