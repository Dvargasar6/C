#include<stdio.h>  //Librería estándar de entrada y salida.
int main(void)
{
    int numero;
    char texto[10]; //Cadena de caracteres de máximo 10 caracteres.
    printf("Ingrese un entero por teclado:\n");
    scanf("%d",&numero);
    printf("Ingrese un texto de máximo 10 caracteres por teclado:\n");
    scanf("%s",texto);
    printf("El numero ingresado es: %d\nEl texto ingresado es: %s\n",numero,texto);
    return(0);

}