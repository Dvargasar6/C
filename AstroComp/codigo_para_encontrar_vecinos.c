#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<malloc.h>

double calcular_distancia_entre_puntos(double x1, double x2, double y1, double y2);
int buscarEnArreglo(const int arreglo[], int busqueda, int longitud);



int main()
{
    FILE *fpuntos, *fvecinos, *fpuntointeres;
    int i,j,k,l; //Contadores.
    int punto_de_interes = 745;  //<-- Puede cambiar el punto de interés aquí.
    int numero_de_puntos = 1000, numero_de_puntos_celdas = 100;
    int dimension=100, numero_de_celdas = dimension*dimension; //Dimensiones de la celda.
    double h = 2*sqrt(2);
    double coordenada_x,coordenada_y; 
    double distancia;
    double puntos[numero_de_puntos][2];
    int celda_del_punto[numero_de_puntos]; 
    double centros_celdas[numero_de_celdas][2];
    int celdas_vecinas[numero_de_celdas][8];
    int puntos_por_celda[numero_de_celdas][numero_de_puntos_celdas];
    int puntos_vecinos[numero_de_puntos][numero_de_puntos_celdas];





    fpuntos = fopen("puntos","w");
    for(i=0;i<numero_de_puntos;i++)
    {
        coordenada_x = drand48()*dimension*1.0;
        coordenada_y = drand48()*dimension*1.0;
        puntos[i][0] = coordenada_x;
        puntos[i][1] = coordenada_y;
        celda_del_punto[i] = dimension*floor(coordenada_y)+ceil(coordenada_x)-1;
        //printf("[%lf,%lf] y su celda es:%d\n",puntos[i][0],puntos[i][1],celda_del_punto[i]);
        fprintf(fpuntos,"%lf %lf\n",coordenada_x,coordenada_y);
    }
    fclose(fpuntos);


    for (i=0; i<numero_de_celdas; i++) 
    {
        coordenada_x = (i % dimension) + 0.5;
        coordenada_y = floor(i / dimension) + 0.5;
        centros_celdas[i][0] = coordenada_x;
        centros_celdas[i][1] = coordenada_y;
        //printf("El centro de la celda %d es: [%lf,%lf]\n",i+1,centros_celdas[i][0],centros_celdas[i][1]);
    }



    for(i=0; i<numero_de_celdas; i++)
    {
        k = 0;
        for(j=0;j<8;j++)
        {
            celdas_vecinas[i][j] = -1;
        }
        for (j=0; j<numero_de_celdas; j++)
        {
            if(i!=j)
            {
                distancia = calcular_distancia_entre_puntos(centros_celdas[i][0],centros_celdas[j][0],
                                                            centros_celdas[i][1],centros_celdas[j][1]);
                if(distancia <= sqrt(2) + 1e-6)  
                {
                    celdas_vecinas[i][k] = j;
                    k++;
                }

            } 
        }
        //printf("Las vecinas de la celda %d son [%d, %d, %d, %d, %d, %d, %d, %d]\n",i,celdas_vecinas[i][0],celdas_vecinas[i][1],celdas_vecinas[i][2],celdas_vecinas[i][3],celdas_vecinas[i][4],celdas_vecinas[i][5],celdas_vecinas[i][6],celdas_vecinas[i][7]);
    }
    for (i=0;i<numero_de_celdas;i++)
    {
        k=0;
        for (j=0;j<numero_de_puntos_celdas;j++)
        {
            puntos_por_celda[i][j]=-1;
        }
        for(j=0;j<numero_de_puntos;j++)
        { 
            //printf("La celda del punto %d es %d\n",j+1,celda_del_punto[j]);
            if((celda_del_punto[j]-1)==i)
            {
                puntos_por_celda[i][k]=j;
                k++;
            }
        }
        //printf("La celda %d tiene los puntos:[%d,%d,%d]\n",i+1,puntos_por_celda[i][0]+1,puntos_por_celda[i][1]+1,puntos_por_celda[i][2]+1);
    }
    for(i=0;i<numero_de_puntos;i++)
    {
        l=0;
        for(j=0;j<numero_de_puntos_celdas;j++)
        {
            puntos_vecinos[i][j]=-1;
        }
        for(j=0;j<numero_de_puntos;j++)
        {
            if (i!=j)
            {
                for(k=0;k<8;k++)
                {   
                    if (calcular_distancia_entre_puntos(puntos[i][0],puntos[j][0],puntos[i][1],puntos[j][1])<=h+1e-10)
                    {    
                        //printf("Celdas vecinas de i:%d, %d\n",i+1,celdas_vecinas[celda_del_punto[i]][k]);
                        if(celda_del_punto[j]==celdas_vecinas[celda_del_punto[i]][k] )
                        {
                            //printf("celda de j:%d, %d celda vecina k: %d\n",j+1,celda_del_punto[j],celdas_vecinas[celda_del_punto[i]][k]);
                            puntos_vecinos[i][l++]=j;
                        }
                        if(celda_del_punto[j]==celda_del_punto[i])
                        {
                            puntos_vecinos[i][l++]=j;
                            break;
                        }
                    }

                }

            }
            
        }
        //printf("Los vecinos de %d son: [%d,%d,%d,%d,%d,%d]\n",i,puntos_vecinos[i][0],puntos_vecinos[i][1],puntos_vecinos[i][2],puntos_vecinos[i][3],puntos_vecinos[i][4],puntos_vecinos[i][5]);
    }
    fvecinos = fopen("archivo_vecinos","w");
    

    for(j=0;j<numero_de_puntos_celdas;j++)
        {
            if(puntos_vecinos[punto_de_interes][j]!=-1)
            {
                //printf("%d\n",puntos_vecinos[9][j]);
                //printf("%lf,%lf\n",puntos[puntos_vecinos[9][j]][0],puntos[puntos_vecinos[9][j]][1]);
                coordenada_x=puntos[puntos_vecinos[punto_de_interes][j]][0];
                coordenada_y=puntos[puntos_vecinos[punto_de_interes][j]][1];
                fprintf(fvecinos,"%lf %lf\n",coordenada_x,coordenada_y);
            }
        }
    fclose(fvecinos);
    fpuntointeres = fopen("punto","w");
    fprintf(fpuntointeres,"%lf %lf\n",puntos[punto_de_interes][0],puntos[punto_de_interes][1]);
    fclose(fpuntointeres);


}


double calcular_distancia_entre_puntos(double x1, double x2, double y1, double y2)
{
    double distancia;
    distancia = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
    return distancia;
}

//Comando de graficación:
//plot "puntos" u 1:2 pt 7 lc rgb "blue" not, "archivo_vecinos" u 1:2 pt 7 lc rgb "orange" not, "punto" u 1:2 pt 7 lc rgb "red" not
