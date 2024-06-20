#pragma once

#ifndef EXAMPLES_REGISTER_H
#define EXAMPLES_REGISTER_H

typedef void (*example_func)();

typedef struct
{
    const char* name;
    example_func func;
} example_entry;

#define REGISTER_EXAMPLE(func) {#func, func}

#endif // EXAMPLES_REGISTER_H
