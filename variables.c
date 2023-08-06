#include<stdio.h>
//#include<stdlib.h>
int main()
{
    float numero_real;   //Creación de la variable.
    numero_real = 34.5;  //Asignación.
    printf("Real: %f\n",numero_real);

    int entero;
    entero = 5;
    printf("Entero: %d\n",entero);

    size_t tamano;
    tamano = sizeof(entero);
    printf("Tamaño de un entero: %lu\n",tamano);

    exit(0); //Salida del programa.

    return 0;
}