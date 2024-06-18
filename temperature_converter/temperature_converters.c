#include "temperature_converters.h"

float celsiusToFahrenheit(float celsius)
{
    return (celsius * 9.0 / 5.0) + 32;
}

float fahrenheitToCelsius(float fahrenheit)
{
    return (fahrenheit - 32) * 5.0 / 9.0;
}