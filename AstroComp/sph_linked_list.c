#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<malloc.h>

#define X 0
#define Y 1

typedef struct
{
  int id;
  double pos[2];
  double vel[2];
  double accel[2];
  double mass;
  double rho;
  double h;
  double p;
  double c;
  double du;
  double u;
  int *nn;
  int nNeighbors;
  double *dx;
  double *dy;
  double *r;
  double *W;
  double *dWx;
  double *dWy;
  int type;
}Particles;

Particles *part, *auxPart;

int nFluid, nPart;

double calcular_distancia_entre_puntos(double x1, double x2, double y1, double y2);
void ics(int nx, int ny, double dx, double dy, double Lx, double Ly);
double W(double r, double h);
double dW(double r, double dx, double h);
void testKernel(void);
void NN(int i);
void test_NN(void);
void density(void);
void eos(void);
void navierStokes(void);
void viscosity(double dx);
void boundaryInteraction(double dx);
void meanVelocity(void);
void acceleration(double dx);
void drift(double dt);
void kick(double dt);
void printState(char *outfile);

int numero_de_puntos = 1600, dimension = 100, numero_de_celdas = 100*100, numero_de_puntos_celdas = 100;
double h = 2*sqrt(2);
int celda_del_punto[1680]; 
int celdas_vecinas[100*100][8];
int puntos_vecinos[1680][100];
int arreglo[8];

int main(int argc, char *argv[])
{

  int i, nx, ny, counter;

  ///////////////////////////////////////////////////
    int w,j,k,l; //Contadores.
    //int punto_de_interes = i;  //<-- Puede cambiar el punto de interés aquí.
    int numero_de_puntos = 1700, numero_de_puntos_celdas = 100;
    int dimension=100, numero_de_celdas = dimension*dimension; //Dimensiones de la celda.
    double h = 2*sqrt(2);
    double coordenada_x,coordenada_y; 
    double distancia;
    double puntos[numero_de_puntos][2];
    double centros_celdas[numero_de_celdas][2];
    int puntos_por_celda[numero_de_celdas][numero_de_puntos_celdas];
    int celda_del_punto[numero_de_puntos]; 
    int celdas_vecinas[numero_de_celdas][8];
    int puntos_vecinos[numero_de_puntos][numero_de_puntos_celdas];
    

  //////////////////////////////////////////77
  
  double Lx, Ly, dx, dy;
  double dt = 5e-5;
  double t, tTotal = atoi(argv[1])*dt;
  char outfiles[500];
  //double t, tTotal = 4000*dt;

  printf("voy a correr durante %d pasos, un tiempo total de %lf s\n",atoi(argv[1]),tTotal);
    
  nx = 40;
  ny = 40;
  Lx = 1e-3;
  Ly = 1e-3;

  dx = Lx/nx;
  dy = Ly/ny;
  
  nFluid = nx*ny;
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
  part = (Particles *)malloc((size_t)nFluid*sizeof(Particles));
  if( part==NULL )
    {
      printf("Error alocando part\n");
      exit(0);
    }
  
  // Create the initial conditions
  ics( nx, ny, dx, dy, Lx, Ly);
 
  // testing kernel function
  testKernel();
       
  counter = 0;
  t = 0;
  
  // printting system initial state
  sprintf(outfiles,"./output/state_%.4d",counter);
  printState(outfiles);




  /////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

      //Este for crea los centros de las celdas

    for (i=0; i<numero_de_celdas; i++) 
    {
        coordenada_x = (i % dimension) + 0.5;
        coordenada_y = floor(i / dimension) + 0.5;
        centros_celdas[i][0] = coordenada_x;
        centros_celdas[i][1] = coordenada_y;
        //printf("El centro de la celda %d es: [%lf,%lf]\n",i+1,centros_celdas[i][0],centros_celdas[i][1]);
    }

    //Este for inicializa las celdas en -1.
    for (w=0;w<numero_de_celdas;w++)
    {
        k=0;
        for (j=0;j<numero_de_puntos_celdas;j++)
        {
            puntos_por_celda[w][j]=-1;
        }
    }

    //Este for crea las celdas vecinas a cada celda
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



    ////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////77
    ////////////////////////////////////////////////////////////////////
  
  // main loop
  while( t<=tTotal )
    {
      
      // searching near neighbors for all fuid particles
      for( i=0; i<nFluid; i++ )
      {

        coordenada_x = part[i].pos[X];
        coordenada_y = part[i].pos[Y];
        //puntos[w][0] = coordenada_x;
        //puntos[w][1] = coordenada_y;
        celda_del_punto[i] = dimension*floor(coordenada_y)+ceil(coordenada_x)-1;
        //printf("[%lf,%lf] y su celda es:%d\n",puntos[i][0],puntos[i][1],celda_del_punto[i]);
        //fprintf(fpuntos,"%lf %lf\n",coordenada_x,coordenada_y);

         // Este for añade los puntos a una celda
        for (k=0;k<numero_de_puntos_celdas;k++)
        {
            if(puntos_por_celda[celda_del_punto[i]-1][k]==-1)
            {
                puntos_por_celda[celda_del_punto[i]-1][k]=i;
            }
            else
            {
                continue;
            }
        }
        for(j=0;j<numero_de_puntos_celdas;j++)
        {
            puntos_vecinos[i][j]=-1;
        }
        NN(i);
      }
      
      // testing near neighbors searching
      if(counter==0)
	test_NN();
      
      // computing density

      density();
            
      // drift in leap-frog integration
      drift(dt);

      // computing acceleration
      acceleration(dx);  

      // kick in leap-frog integration
      kick(dt);
      
      // drift in leap-frog integration
      drift(dt);
	
      t = t + dt;
      counter++;

      // printting system state
      sprintf(outfiles,"./output/state_%.4d",counter);
      printState(outfiles);

      printf("step = %d \n",counter);
      
    }
  
  free(part);
  
  return 0;
}

void ics(int nx, int ny, double dx, double dy, double Lx, double Ly)
{
  int i, j, counter;
  
  FILE *fFluidIcs, *fbBorder, *frBorder, *ftBorder, *flBorder;
  fFluidIcs = fopen("fluid_ics.output","w");
    
  // ics for fluid particles
  
  counter = 0;
  for( j=0; j<ny; j++)
    {
      for( i=0; i<nx; i++)
	{
	  part[counter].id = counter;
	  part[counter].pos[X] = i*dx+dx/2.0;
	  part[counter].pos[Y] = j*dy+dy/2.0;
	  part[counter].vel[X] = 0.0;
	  part[counter].vel[Y] = 0.0;
	  part[counter].accel[X] = 0.0;
	  part[counter].accel[Y] = 0.0;
	  part[counter].rho = 1000;
	  part[counter].h = dx;
	  part[counter].mass = part[counter].rho*dx*dy;
	  part[counter].p = 0.0;
	  part[counter].c = 0.0;
	  part[counter].du = 0.0;
	  part[counter].u = 357.1; 
	  part[counter].nn = NULL;
	  part[counter].nNeighbors = 0;
	  part[counter].dx = NULL;
	  part[counter].dy = NULL;
	  part[counter].r = NULL;
	  part[counter].W = NULL;
	  part[counter].dWx = NULL;
	  part[counter].dWy = NULL;
	  part[counter].type = 1;
	  counter++;
	}
    }
  
  // ics for boundary particles

  // speed in boundary
  double vBoundary = 1.5e-2; 
  
  int npVirtI = 320;
  int npV = npVirtI/4;
  
  // bottom border, 81 points

  fbBorder = fopen("bottom_border.output","w");

  //Boundary particles are add to estructure particles set
  nPart = nFluid;
  
  auxPart = NULL;
  
  auxPart = (Particles *)realloc(part, (size_t)(nPart+npV+1)*sizeof(Particles));
  if(auxPart==NULL)
    {
      printf("error en auxPart\n");
      exit(0);
    }
  else
    {
      part = auxPart;
      auxPart = NULL;
    }
  
  counter = nPart;
    
  for( i=0; i<=npV; i++)
    {
      part[counter].id = counter;
      part[counter].pos[X] = i*dx/2.0;
      part[counter].pos[Y] = 0.0;
      part[counter].vel[X] = 0.0;
      //part[counter].vel[X] = -vBoundary;
      part[counter].vel[Y] = 0.0;
      part[counter].accel[X] = 0.0;
      part[counter].accel[Y] = 0.0;
      part[counter].rho = 1000;
      part[counter].h = dx;
      part[counter].mass = part[counter].rho*dx*dy;
      part[counter].p = 0.0;
      part[counter].c = 0.0;
      part[counter].du = 0.0;
      part[counter].u = 357.1;
      part[counter].nn = NULL;
      part[counter].nNeighbors = 0;
      part[counter].dx = NULL;
      part[counter].dy = NULL;
      part[counter].r = NULL;
      part[counter].W = NULL;
      part[counter].dWx = NULL;
      part[counter].dWy = NULL;
      part[counter].type = -1;
      counter++;
    }

   for( i=nPart; i<nPart+npV+1; i++)
    {
      fprintf(fbBorder,"%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
	      part[i].id,
	      part[i].pos[X],part[i].pos[Y],
	      part[i].vel[X],part[i].vel[Y],
	      part[i].accel[X],part[i].accel[Y],
	      part[i].rho,part[i].mass,
	      part[i].p,part[i].c,part[i].u);
    }
  
  fclose(fbBorder);
  
  // right border, 79 points

  frBorder = fopen("right_border.output","w");

  //Boundary particles are add to estructure particles set
  nPart = counter;
  
  auxPart = NULL;
    
  auxPart = (Particles *)realloc(part,(size_t)(nPart+npV-1)*sizeof(Particles));
  if(auxPart==NULL)
    {
      printf("error en auxPart\n");
      exit(0);
    }
  else
    {
      part = auxPart;
      auxPart = NULL;
    }
    
  for( i=0; i<npV-1; i++)
    {
      part[counter].id = counter;
      part[counter].pos[X] = Lx;
      part[counter].pos[Y] = dy/2.0 + i*dy/2.0;
      part[counter].vel[X] = 0.0;
      part[counter].vel[Y] = 0.0;
      //part[counter].vel[Y] = -vBoundary;
      part[counter].accel[X] = 0.0;
      part[counter].accel[Y] = 0.0;
      part[counter].rho = 1000;
      part[counter].h = dx;
      part[counter].mass = part[counter].rho*dx*dy;
      part[counter].p = 0.0;
      part[counter].c = 0.0;
      part[counter].du = 0.0;
      part[counter].u = 357.1;
      part[counter].nn = NULL;
      part[counter].nNeighbors = 0;
      part[counter].dx = NULL;
      part[counter].dy = NULL;
      part[counter].r = NULL;
      part[counter].W = NULL;
      part[counter].dWx = NULL;
      part[counter].dWy = NULL;
      part[counter].type = -1;
      counter++;
    }
  
    for( i=nPart; i<nPart+npV-1; i++)
    {
      fprintf(frBorder,"%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
	      part[i].id,
	      part[i].pos[X],part[i].pos[Y],
	      part[i].vel[X],part[i].vel[Y],
	      part[i].accel[X],part[i].accel[Y],
	      part[i].rho,part[i].mass,
	      part[i].p,part[i].c,part[i].u);
    }
  
  fclose(frBorder);
  
  // top border, 81 points

  ftBorder = fopen("top_border.output","w");


  //Boundary particles are add to estructure particles set
  nPart = counter;
  
  auxPart = NULL;
  
  auxPart = (Particles *)realloc(part,(size_t)(nPart+npV+1)*sizeof(Particles));
  if(auxPart==NULL)
    {
      printf("error en auxPart\n");
      exit(0);
    }
  else
    {
      part = auxPart;
      auxPart = NULL;
    }
      
  for( i=0; i<=npV; i++)
    {
      part[counter].id = counter;
      part[counter].pos[X] = i*dx/2.0;
      part[counter].pos[Y] = Ly;
      part[counter].vel[X] = vBoundary;
      part[counter].vel[Y] = 0.0;
      part[counter].accel[X] = 0.0;
      part[counter].accel[Y] = 0.0;
      part[counter].rho = 1000;
      part[counter].h = dx;
      part[counter].mass = part[counter].rho*dx*dy;
      part[counter].p = 0.0;
      part[counter].c = 0.0;
      part[counter].du = 0.0;
      part[counter].u = 357.1;
      part[counter].nn = NULL;
      part[counter].nNeighbors = 0;
      part[counter].dx = NULL;
      part[counter].dy = NULL;
      part[counter].r = NULL;
      part[counter].W = NULL;
      part[counter].dWx = NULL;
      part[counter].dWy = NULL;
      part[counter].type = -1;
      counter++;
    }

  for( i=nPart; i<nPart+npV+1; i++)
    {
      fprintf(ftBorder,"%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
	      part[i].id,
	      part[i].pos[X],part[i].pos[Y],
	      part[i].vel[X],part[i].vel[Y],
	      part[i].accel[X],part[i].accel[Y],
	      part[i].rho,part[i].mass,
	      part[i].p,part[i].c,part[i].u);
    }
  
  fclose(ftBorder);
  
  // left border, 79 points

  flBorder = fopen("left_border.output","w");
  

  //Boundary particles are add to estructure particles set
  nPart = counter;
  
  auxPart = NULL;
    
  auxPart = (Particles *)realloc(part,(size_t)(nPart+npV-1)*sizeof(Particles));
  if(auxPart==NULL)
    {
      printf("error en auxPart\n");
      exit(0);
    }
  else
    {
      part = auxPart;
      auxPart = NULL;
    }
  
  for( i=0; i<npV-1; i++)
    {
      part[counter].id = counter;
      part[counter].pos[X] = 0.0;
      part[counter].pos[Y] = dy/2.0 + i*dy/2.0;
      part[counter].vel[X] = 0.0;
      part[counter].vel[Y] = 0.0;
      // part[counter].vel[Y] = vBoundary;
      part[counter].accel[X] = 0.0;
      part[counter].accel[Y] = 0.0;
      part[counter].rho = 1000;
      part[counter].h = dx;
      part[counter].mass = part[counter].rho*dx*dy;
      part[counter].p = 0.0;
      part[counter].c = 0.0;
      part[counter].du = 0.0;
      part[counter].u = 357.1;
      part[counter].nn = NULL;
      part[counter].nNeighbors = 0;
      part[counter].dx = NULL;
      part[counter].dy = NULL;
      part[counter].r = NULL;
      part[counter].W = NULL;
      part[counter].dWx = NULL;
      part[counter].dWy = NULL;
      part[counter].type = -1;
      counter++;
    }

   for( i=nPart; i<nPart+npV-1; i++)
    {
      fprintf(flBorder,"%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
	      part[i].id,
	      part[i].pos[X],part[i].pos[Y],
	      part[i].vel[X],part[i].vel[Y],
	      part[i].accel[X],part[i].accel[Y],
	      part[i].rho,part[i].mass,
	      part[i].p,part[i].c,part[i].u);
    }
  
  fclose(flBorder);
  
  // print all particles

  nPart = counter;
 
  for( i=0; i<nPart; i++)
    {
      fprintf(fFluidIcs,"%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
	      part[i].id,
	      part[i].pos[X],part[i].pos[Y],
	      part[i].vel[X],part[i].vel[Y],
	      part[i].accel[X],part[i].accel[Y],
	      part[i].rho,part[i].mass,
	      part[i].p,part[i].c,part[i].u);
    }

  fclose(fFluidIcs);

  
}

double W(double r, double h)
{
  
  double R = r/h;
  
  double alpha = 15.0/(7.0*M_PI*h*h);
  
  if( (R >= 0.0) && (R < 1.0) )
    return alpha*((2.0/3.0) - R*R + 0.5*R*R*R);
  
  if( (R >= 1.0) && (R <= 2.0) )
    return alpha*((1.0/6.0)*(2.0-R)*(2.0-R)*(2.0-R));

  if( R>2.0)
    return 0.0;
  
  return 0.0;
}

double dW(double r, double dx, double h)
{
  
  double R = r/h;
  
  double alpha = 15.0/(7.0*M_PI*h*h);
  
  if( (R >= 0.0) && (R < 1.0) )
    return alpha*(-2.0 + 1.5*R)*dx/(h*h);
  
  if( (R >= 1.0) && (R <= 2.0) )
    return alpha*(-0.5*(2.0-R)*(2.0-R))*dx/(h*h*R);

  if( R>2.0)
    return 0.0;
  
  return 0.0;
}

void testKernel(void)
{
  double r, w, dw;

  FILE *fKernelTest;
  fKernelTest = fopen("kernel_test.output","w");
  
  for( r=-3.0; r<=3.0; r = r + 0.1)
    {
      w = W( fabs(r), 1.0);
      dw = dW( fabs(r), r/sqrt(3.0), 1.0);

      fprintf(fKernelTest,"%16.10lf %16.10lf %16.10lf\n",r,w,dw);
      
    }

  fclose(fKernelTest);
  
}

// Searching the near neighbors 
void NN(int i)
{

    double kappa = 2.0;
  double xij, yij, rij, hij;
  double *auxDouble;
  int j, *auxInt, nNeighbors, k,l;
   
  free(part[i].nn);
  free(part[i].dx);
  free(part[i].dy);
  free(part[i].r);
  free(part[i].W);
  free(part[i].dWx);
  free(part[i].dWy);

  part[i].nn = NULL;
  part[i].dx = NULL;
  part[i].dy = NULL;
  part[i].r = NULL;
  part[i].W = NULL;
  part[i].dWx = NULL;
  part[i].dWy = NULL;



  ///////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////



    //Este for añade los vecinos de cada punto.


        l=0;
        nNeighbors = 0;
        for(j=0;j<numero_de_puntos;j++)
        hij = 0.5*(part[i].h+part[j].h);
        {
            if (i!=j)
            {
                for(k=0;k<8;k++)
                {   
                    if (calcular_distancia_entre_puntos(part[i].pos[X],part[j].pos[X],part[i].pos[Y],part[j].pos[Y])<=kappa*hij+1e-10)
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

        nNeighbors++;
                        	      // add neighbor id
	      auxInt = NULL;
	      auxInt = (int *)realloc(part[i].nn,(size_t)(nNeighbors)*sizeof(int));
	      part[i].nn = auxInt;
	      auxInt = NULL;
	      	      
	      // add neighbor dx
	      auxDouble = NULL;
	      auxDouble = (double *)realloc(part[i].dx,(size_t)(nNeighbors)*sizeof(double));
	      part[i].dx = auxDouble;
	      auxDouble = NULL;

	      // add neighbor dy
	      auxDouble = NULL;
	      auxDouble = (double *)realloc(part[i].dy,(size_t)(nNeighbors)*sizeof(double));
	      part[i].dy = auxDouble;
	      auxDouble = NULL;

	      // add neighbor r
	      auxDouble = NULL;
	      auxDouble = (double *)realloc(part[i].r,(size_t)(nNeighbors)*sizeof(double));
	      part[i].r = auxDouble;
	      auxDouble = NULL;
	      
	      // add neighbor W
	      auxDouble = NULL;
	      auxDouble = (double *)realloc(part[i].W,(size_t)(nNeighbors)*sizeof(double));
	      part[i].W = auxDouble;
	      auxDouble = NULL;

	      // add neighbor dWx
	      auxDouble = NULL;
	      auxDouble = (double *)realloc(part[i].dWx,(size_t)(nNeighbors)*sizeof(double));
	      part[i].dWx = auxDouble;
	      auxDouble = NULL;
	      
	      // add neighbor dWy
	      auxDouble = NULL;
	      auxDouble = (double *)realloc(part[i].dWy,(size_t)(nNeighbors)*sizeof(double));
	      part[i].dWy = auxDouble;
	      auxDouble = NULL;
	      
	      part[i].nn[nNeighbors-1] = j;
	      part[i].dx[nNeighbors-1] = xij;
	      part[i].dy[nNeighbors-1] = yij;
	      part[i].r[nNeighbors-1] = rij;
	      part[i].W[nNeighbors-1] = W( rij, hij ); 
	      part[i].dWx[nNeighbors-1] = dW( rij, xij, hij);
	      part[i].dWy[nNeighbors-1] = dW( rij, yij, hij);
                    }

                }

            }
            
        }
        //printf("Los vecinos de %d son: [%d,%d,%d,%d,%d,%d]\n",i,puntos_vecinos[i][0],puntos_vecinos[i][1],puntos_vecinos[i][2],puntos_vecinos[i][3],puntos_vecinos[i][4],puntos_vecinos[i][5]);


    ///////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////

  part[i].nNeighbors = nNeighbors;











}

double calcular_distancia_entre_puntos(double x1, double x2, double y1, double y2)
{
    double distancia;
    distancia = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
    return distancia;
}

void test_NN(void)
{
  int i,j,k;

   FILE *fTestNN;
   fTestNN = fopen("NN_test.output","w");
   srand(time(NULL));
   
   for( k=0; k<20; k++)
    {

      i = rand() % nFluid;
      
      printf("testing for particle %d\n",i);
      printf("with %d neighbors\n",part[i].nNeighbors);
       
      fprintf(fTestNN,"%16d %16.10lf %16.10lf\n",
	      part[i].id,
	      part[i].pos[X],
	      part[i].pos[Y]);
      
      for( j=0; j<part[i].nNeighbors; j++ )
	fprintf(fTestNN,"%16d %16.10lf %16.10lf\n",
		part[i].nn[j],
		part[part[i].nn[j]].pos[X],
		part[part[i].nn[j]].pos[Y]);
      fprintf(fTestNN,"\n");
    }
  fclose(fTestNN);
    
}

void density(void)
{
  int i, j;
  double wii, norm;
  
  for( i=0; i<nFluid; i++ )
    {
      // self density
      wii = W( 0.0, part[i].h );
      
      // computing density
      part[i].rho = part[i].mass*wii;
      for( j=0; j<part[i].nNeighbors; j++ )
	part[i].rho = part[i].rho + part[part[i].nn[j]].mass*part[i].W[j];

      // normalizing the density
      norm = (part[i].mass/part[i].rho)*wii;
      for( j=0; j<part[i].nNeighbors; j++ )
	  norm = norm + (part[part[i].nn[j]].mass/part[part[i].nn[j]].rho)*part[i].W[j];
      
      part[i].rho = part[i].rho/norm;
    }
  
  printf("density computed\n");

}

void eos(void)
{
  int i;
  
  for( i=0; i<nPart; i++ )
    {
      part[i].c = 0.01;
      part[i].p = part[i].c*part[i].c*part[i].rho; 
    }
}

void navierStokes(void)
{

  int i, j, k;
  double pij, vdw;
  // computing sound speed and pression
  eos();

  // computing acceleration
  for( i=0; i<nFluid; i++ )
    {
      
      part[i].accel[X] = part[i].accel[Y] = 0.0;
      part[i].du = 0.0;

      for( k=0; k<part[i].nNeighbors; k++ )
	{
	  j = part[i].nn[k];
	  
	  pij = ( part[i].p/(part[i].rho*part[i].rho) )
	    + ( part[j].p/(part[j].rho*part[j].rho) );
	  part[i].accel[X] = part[i].accel[X] - part[j].mass*pij*part[i].dWx[k];
	  part[i].accel[Y] = part[i].accel[Y] - part[j].mass*pij*part[i].dWy[k];
	  
	  vdw = (part[i].vel[X]-part[j].vel[X])*part[i].dWx[k]
	    + (part[i].vel[Y]-part[j].vel[Y])*part[i].dWy[k];
	  part[i].du = part[i].du + 0.5*part[j].mass*pij*vdw;
	  
	} //W = integral( f*dr ) = integral( f*dr*(dt/dt) ) = integral( f*v*dt )

    }
  
  printf("acceleration computed\n");
}

void viscosity(double dx)
{
  
  int i, j, k;

  double xij, yij, vxij, vyij, vijrij, vdw;
  double hij, cij, phiij, rhoij, Piij;
  double alphapi = 1.0;
  double betapi = 1.0;
  double eps = dx;
  double eps2 = 0.01*eps*eps;

  for( i=0; i<nFluid; i++ )
    {
      for( k=0; k<part[i].nNeighbors ; k++ )
	{
	  
	  j = part[i].nn[k];
	  
	  xij = part[i].pos[X] - part[j].pos[X];
	  yij = part[i].pos[Y] - part[j].pos[Y];
	  vxij = part[i].vel[X] - part[j].vel[X];
	  vyij = part[i].vel[Y] - part[j].vel[Y];
	  vijrij = vxij*xij + vyij*yij;
	  
	  if( vijrij < 0.0 )
	    {
	      hij = 0.5*(part[i].h+part[j].h);
	      phiij = (hij*vijrij)/( xij*xij + yij*yij + eps2);
	      cij = 0.5*(part[i].c+part[j].c);
	      rhoij = 0.5*(part[i].rho+part[j].rho);
	      
	      Piij = ( -alphapi*cij*phiij + betapi*phiij*phiij )/( rhoij );

	      part[i].accel[X] = part[i].accel[X] - part[j].mass*Piij*part[i].dWx[k];
	      part[i].accel[Y] = part[i].accel[Y] - part[j].mass*Piij*part[i].dWy[k];

	      vdw = (part[i].vel[X]-part[j].vel[X])*part[i].dWx[k]
		+ (part[i].vel[Y]-part[j].vel[Y])*part[i].dWy[k];
	      part[i].du = part[i].du + 0.5*part[j].mass*Piij*vdw;
	      
	    }
	  
	}
    }
  printf("viscosity computed\n");
}

void boundaryInteraction(double dx)
{

  int i, j;
  int n1 = 12, n2 = 4;
  double r0 = dx/2.0, D = 0.01;
  double xij, yij, rij, PBxij, PByij;

  for( i=0; i<nFluid; i++ )
    {
      for( j=0; j<part[i].nNeighbors; j++ )
	{
	  if( part[part[i].nn[j]].type==-1 )
	    {
	      xij = part[i].pos[X] - part[part[i].nn[j]].pos[X];
	      yij = part[i].pos[Y] - part[part[i].nn[j]].pos[Y];
	      rij = sqrt( xij*xij + yij*yij );
	      
	      if( rij<r0 )
		{
		  PBxij = D*( pow((r0/rij),n1) - pow((r0/rij),n2) )*(xij/(rij*rij));
		  PByij = D*( pow((r0/rij),n1) - pow((r0/rij),n2) )*(yij/(rij*rij));
		  
		  part[i].accel[X] = part[i].accel[X] + PBxij;
		  part[i].accel[Y] = part[i].accel[Y] + PByij;
		  
		}
	    }
	}
    }
  printf("interaction with boundary computed\n");
}

void meanVelocity(void)
{

  int i, j;
  double epsilon = 0.3;
  double vxMean, vyMean;
  double vxij, vyij, rhoij;

  for( i=0; i<nFluid; i++ )
    {
      
      vxMean = 0.0;
      vyMean = 0.0;
      
      for( j=0; j<part[i].nNeighbors; j++ )
	{
	  vxij = part[i].vel[X] - part[part[i].nn[j]].vel[X];
	  vyij = part[i].vel[Y] - part[part[i].nn[j]].vel[Y];
	  rhoij = 0.5*(part[i].rho+part[part[i].nn[j]].rho);
	  vxMean = vxMean + (part[part[i].nn[j]].mass/rhoij)*vxij*part[i].W[j];
	  vyMean = vyMean + (part[part[i].nn[j]].mass/rhoij)*vyij*part[i].W[j];
	}

      part[i].vel[X] = part[i].vel[X] - epsilon*vxMean;
      part[i].vel[Y] = part[i].vel[Y] - epsilon*vyMean;
   
    }
}

void acceleration(double dx)
{

  // computing acceleration and change of energy
  navierStokes();
  
  // computing viscosity contribution
  viscosity(dx);
  
  // computing interaction with boundary
  boundaryInteraction(dx);

  // correction to mean velocity
  meanVelocity();

  printf("acceleration computed\n");
  
}

void drift(double dt)
{

  int i;
  
  for( i=0; i<nFluid; i++ )
    {
      part[i].pos[X] = part[i].pos[X] + 0.5*dt*part[i].vel[X];
      part[i].pos[Y] = part[i].pos[Y] + 0.5*dt*part[i].vel[Y];
      part[i].u = part[i].u + 0.5*dt*part[i].du;
    }
}

void kick(double dt)
{
  
  int i;
  
  for( i=0; i<nFluid; i++ )
    {
      part[i].vel[X] = part[i].vel[X] + dt*part[i].accel[X];
      part[i].vel[Y] = part[i].vel[Y] + dt*part[i].accel[Y];
    }
}

void printState(char *outfile)
{

  int i;

  FILE *fState;
  fState = fopen(outfile,"w");
 
  for( i=0; i<nPart; i++)
    {
      fprintf(fState,"%d %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
	      part[i].id,
	      part[i].pos[X],part[i].pos[Y],
	      part[i].vel[X],part[i].vel[Y],
	      part[i].accel[X],part[i].accel[Y],
	      part[i].rho,part[i].mass,
	      part[i].p,part[i].c,part[i].u);
    }
  
  fclose(fState);

}
