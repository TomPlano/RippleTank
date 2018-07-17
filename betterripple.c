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

if(argc!=5)
{
    printf("%i\n",argc);
    printf("./ripple.x <m> <n> <time> <itters>\n");
    return 0;
}


int M= atoi(argv[1]);
int N= atoi(argv[1]);
int T= atoi(argv[3]);
int L= atoi(argv[4]); 

double dx=1.0/(M+1.0);
double dy=1.0/(N+1.0);
double dt=(1.0*T)/(1.0*L);

double** grid_ptr= (double**)allocMat(M,sizeof(double));//sets all to zero
double** next = (double**)allocMat(M,sizeof(double));
double** prior = (double**)allocMat(M,sizeof(double));

double** temp;


/* Enforcing initial condition 1 */
for (int j=0; j<N; j++) {
double y = j*dy;
for (int i=0; i<M; i++) {
double x = i*dx;
prior[j][i] = sin(2*M_PI*(x+y));
}
}

/* Enforcing initial condition 2 */
for (int j=1; j<N-1; j++)
for (int i=1; i<M-1; i++) /* interior points */
grid_ptr[j][i] = prior[j][i]+((dt*dt)/(4*dx*dx))*(prior[j][i-1]-2*prior[j][i]+prior[j][i+1])+((dt*dt)/(4*dy*dy))*(prior[j-1][i]-2*prior[j][i]+prior[j+1][i]);

#pragma omp parallel
{
    for(double timestep=0; timestep<T; timestep+=dt)
    {
            
            #pragma omp for schedule(dynamic) 
            for(int i=1;i<N-1;i++)
            {
                    for(int j=1;j<M-1;j++)
                    {
                        next[i][j]=discrete_wave_eq(grid_ptr,prior,dx,dy,dt,i,j);
                    }
            }
        #pragma omp master
        {
            //print out cuttent ittereation
            pgrid(grid_ptr,M);
            printf("\n");
            //
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
            printf("%f\t",grid[i][j]);
        }
        printf("\n");
    } 
}


void clear(double** grid,int size)
{
  for (int i = 0; i<size-1; i++)
  memset(grid[i],0,size);
  
}
