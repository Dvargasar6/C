#include<stdio.h>
#include<stdlib.h>
int main()
{
    int i; //Es necesario definir la varible de iteración.
    //Ciclo FOR:
    printf("Ciclo for:\n");
    for(i=0;i<=5;i=i+1)
    {
        printf("%d\n",i);
    }

    //Ciclo WHILE:
    printf("Ciclo while:\n");
    while(i>=0)
    {
        printf("%d\n",i);
        i--;
    }
    
    //Ciclo DO...WHILE:
    printf("Ciclo do...while:\n");
    do {printf("%d\n",i); i++;}
    while(i<7);


    return(0);
}