#pragma once

#ifndef EXAMPLES_H
#define EXAMPLES_H

#include "examples_register.h"

// Declare your example functions
void complex_numbers();
void matrix();

// Register example functions
#define EXAMPLES_LIST \
    X(complex_numbers) \
    X(matrix)

#define X(func) REGISTER_EXAMPLE(func),

#define REGISTER_EXAMPLES EXAMPLES_LIST

#endif // EXAMPLES_H