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
//TODO: Change the user input to something better and will work with mpi
if(argc!=4)
{
    printf("./ripple.x <size> <time> <itters>\n");
    return 0;
}

//TODO: support non square sufaces at some point
int M= atoi(argv[1]);//Size of surface, will be a square
int T= atoi(argv[2]);//Total time of simulation
int L= atoi(argv[3]);//Number of itterations, will be used to compute a stable dt eventually

double dx=1.0/(M+1.0);
double dy=1.0/(M+1.0);
//TODO: figure out how to stablize dt for all sizes
double dt=0.015;

//allocMat actually calls calloc, thus is zeroed
double** grid_ptr= (double**)allocMat(M+2,sizeof(double));
double** next = (double**)allocMat(M+2,sizeof(double));
double** prior = (double**)allocMat(M+2,sizeof(double));

double** temp;

#pragma omp parallel
{
    for(double timestep=0; timestep<T; timestep+=dt)
    {
            #pragma omp for
            for(int i=0;i<M+2;i++) //TODO: remove these hardcoded M's at some point
            {
                for(int j=0;j<M+2;j++)
                {
                    //first if statment defines our osilator
                    //TODO: find a better what to define osilators
                    if (i==j&&(i==M/4 || i==3*(M/4)))  {
                        next[i][j] = 5*sin(2*M_PI*3*timestep);
                    }

                    //defines edge behavior, meant to implement Neumann bounday condition
                    //TODO: find better way implementing this, as well as a way to change to different bounds condition
                    else if (i==0){
                        next[i][j] = grid_ptr[i+1][j];
                    }
                    else if (j==0){
                        next[i][j] = grid_ptr[i][j+1];
                    }
                    else if (i==M+1){
                        next[i][j] = grid_ptr[i-1][j];
                    }
                    else if (j==M+1){
                        next[i][j] = grid_ptr[i][j-1];
                    }

                    //main body of work from PDE's for all non-boundry, non-osilators
                    else{
                        next[i][j]=discrete_wave_eq(grid_ptr,prior,dx,dy,dt,i,j);
                    }
                }
            }

        //TODO:
        #pragma omp master
        {
            //print grid in gnuplot readable format
            pgrid(grid_ptr,M+2);
            printf("\n");
            printf("\n");
            //clear old data no longer used
            clear(prior,M);
            //ptr swap for next time step
            temp = next;
            next=prior;
            prior = grid_ptr;
            grid_ptr = temp;
        }
        #pragma omp barrier
    }
}
//cleanup
freeMat(grid_ptr,M);
freeMat(next,M);
freeMat(prior,M);
return 0;
}


double discrete_wave_eq(double** u , double** um, double dx, double dy, double dt, int i, int j){
    //bounds checking is done in loop, this will go out of bounds if you let it
    return  2*u[i][j]-um[i][j]+((dt*dt)/(2*dx*dx))*(u[i][j-1]-2*u[i][j]+u[i][j+1])+((dt*dt)/(2*dy*dy))*(u[i-1][j]-2*u[i][j]+u[i+1][j]);
}


//TODO: support non square sufaces at some point
void pgrid(double** grid, int size){
    for(int i=0;i<size;i++){
        for(int j=0;j<size;j++){
            printf("%.*f ",10 ,grid[i][j]);
        }
        printf("\n");
    } 
}


//TODO: support non square sufaces at some point
void clear(double** grid,int size){
  for (int i = 0; i<size-1; i++)
  memset(grid[i],0,size);
  
}
