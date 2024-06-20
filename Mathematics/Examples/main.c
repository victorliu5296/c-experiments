#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "examples.h"

example_entry examples[] = {
    #define X(func) REGISTER_EXAMPLE(func),
    REGISTER_EXAMPLES
    #undef X
};

void clear_screen()
{
#ifdef _WIN32
    system("cls");
#else
    system("clear");
#endif
}

static void show_menu()
{
    printf("Available examples:\n");
    for (size_t i = 0; i < sizeof(examples) / sizeof(example_entry); ++i)
    {
        printf("%zu: %s\n", i + 1, examples[i].name);
    }
    printf("0: Exit\n");
    printf("-1: Clear Screen\n\n");
}

int main()
{
    while (1)
    {
        show_menu();

        printf("Enter the number of the example to run, -1 to clear the screen, or 0 to exit: ");
        int choice;
        scanf_s("%d", &choice);

        printf("\n");

        if (choice == 0)
        {
            break;
        }
        else if (choice == -1)
        {
            clear_screen();
        }
        else if (choice > 0 && choice <= sizeof(examples) / sizeof(example_entry))
        {
            examples[choice - 1].func();
            printf("\n\n");
        }
        else
        {
            printf("Invalid choice. Please try again.\n\n");
        }
    }

    return 0;
}
