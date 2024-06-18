#include <stdio.h>
#include "temperature_converters.h"

int main()
{
    float temp;
    char unit;

    printf("Enter temperature followed by unit (C/F), e.g. 20 C: ");
    scanf_s("%f %c", &temp, &unit, 1);

    if (unit == 'C' || unit == 'c')
    {
        printf("%.2f C = %.2f F\n", temp, celsiusToFahrenheit(temp));
    }
    else if (unit == 'F' || unit == 'f')
    {
        printf("%.2f F = %.2f C\n", temp, fahrenheitToCelsius(temp));
    }
    else
    {
        printf("Invalid unit.\n");
    }

    return 0;
}