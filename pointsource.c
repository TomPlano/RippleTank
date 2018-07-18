#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <math.h>
#include "Utils.h"

void pgrid(double** grid, int size);
void clear(double** grid,int size);
double discrete_wave_eq(double** u , double** um, double dx, double dy, double dt, int i, int j);
int main(int argc, char** argv)
{

if(argc!=4)
{
    printf("./ripple.x <size> <time> <itters>\n");
    return 0;
}
int M= atoi(argv[1]);
int T= atoi(argv[2]);
int L= atoi(argv[3]); 

double dx=1.0/(M+1.0);
double dy=1.0/(M+1.0);
double dt=0.015;
if(dt>0.015)
{
    printf("please choose <time> and <itters> such that <time>/<itters> < 0.015\n");
    return 0;
}
double** grid_ptr= (double**)allocMat(M+2,sizeof(double));//sets all to zero
double** next = (double**)allocMat(M+2,sizeof(double));
double** prior = (double**)allocMat(M+2,sizeof(double));

double** temp;



#pragma omp parallel
{
    for(double timestep=0; timestep<T; timestep+=dt)
    {
            
            #pragma omp for
            for(int i=0;i<=M+1;i++)
            {
                    for(int j=0;j<=M+1;j++)
                    {
                        if (i==0 || i==M+1 || j==0||j==M+1)
                        {//boundry condition
                            switch(i) {
                                   case 0:
                                      next[i][j]=grid_ptr[i+1][j]; 
                                      break; 
                                   case 51 :
                                      next[i][j]=grid_ptr[i-1][j]; 
                                      break; 
                            }
                            switch(j) {
                                   case 0:
                                      next[i][j]=next[i][j+1]; 
                                      break; 
                                   case 51 :
                                      next[i][j]=next[i][j-1]; 
                                      break; 
                            }
                        }
                        else if (i==j){//line source
                            next[i][j] = 5*sin(2*M_PI*3*timestep);
                        }
                        else{
                            next[i][j]=discrete_wave_eq(grid_ptr,prior,dx,dy,dt,i,j);
                        }
                    }
            }
        #pragma omp master
        {
            pgrid(grid_ptr,M+2);
            printf("\n");
            printf("\n");
            clear(prior,M);
            temp = next;
            next=prior;
            prior = grid_ptr;
            grid_ptr = temp;
    
        }
        #pragma omp barrier
    }
}

freeMat(grid_ptr,M);
freeMat(next,M);
freeMat(prior,M);
return 0;
}

double discrete_wave_eq(double** u , double** um, double dx, double dy, double dt, int i, int j)
{

return  2*u[i][j]-um[i][j]+((dt*dt)/(2*dx*dx))*(u[i][j-1]-2*u[i][j]+u[i][j+1])+((dt*dt)/(2*dy*dy))*(u[i-1][j]-2*u[i][j]+u[i+1][j]);


}



void pgrid(double** grid, int size)
{
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<size;j++)
        {
            printf("%.*f ",10 ,grid[i][j]);

        }
        printf("\n");
    } 
}


void clear(double** grid,int size)
{
  for (int i = 0; i<size-1; i++)
  memset(grid[i],0,size);
  
}
