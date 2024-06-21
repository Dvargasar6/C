#include <stdio.h>

int main(void)
{
    // Variables types:
    int integer;
    float floating;
    double doubleFloat;

    printf("Input a integer:\n");
    scanf("%d", &integer);
    printf("Input a float:\n");
    scanf("%f", &floating);
    printf("Input a double float:\n");
    scanf("%lf", &doubleFloat);

    printf("The integer is: %d\n \
     the float is: %f\n \
      the double float is: %lf\n"\
      ,integer,floating,doubleFloat);

    // Variables size:
    size_t sizeInt, sizeFloat, sizeDouble;
    sizeInt = sizeof(integer);
    sizeFloat = sizeof(floating);
    sizeDouble = sizeof(doubleFloat);

    printf("The size of integer is: %lu\n \
    the size of float is: %lu\n \
    the size of double float is: %lu"\
    ,sizeInt,sizeFloat,sizeDouble);

    return 0;
}
