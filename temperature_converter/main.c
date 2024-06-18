#include <stdio.h>
#include "temperature_converters.h"

int main()
{
    float temperature;
    char unit;

    printf("Enter temperature followed by unit (C/F): ");
    scanf("%f %c", &temperature, &unit);

    switch (unit)
    {
    case 'C':
    case 'c':
        printf("%.2f C = %.2f F\n", temperature, celsiusToFahrenheit(temperature));
        break;
    case 'F':
    case 'f':
        printf("%.2f F = %.2f C\n", temperature, fahrenheitToCelsius(temperature));
        break;
    default:
        printf("Invalid unit.\n");
        break;
    }
    return 0;
}