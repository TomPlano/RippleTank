#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include "Utils.h"

void pgrid(double** grid, int size);
void clear(double** grid,int size);
double discrete_wave_eq(double** u , double** um, double dx, double dy, double dt, int i, int j);
void send_to_printer(double**grid, int size, MPI_Comm comm);

int main(int argc, char** argv)
{
int numranks, rank;
MPI_Status Stat;
MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &numranks);

//TODO: support non square sufaces at some point
int M = 60;//Size of surface, will be a square
int T= 4;//Total time of simulation
int L= 1;//Number of itterations, will be used to compute a stable dt eventually

double dx=1.0/(M+1.0);
double dy=1.0/(M+1.0);
//TODO: figure out how to stablize dt for all sizes
double dt=0.015;


//numranks must eq 17, 0 is printer
int size = M/4;
//taurous world cause reasons
int rank_layout[6][6]={{0,13,14,15,16,0},
                      {4,1,2,3,4,1},
                      {8,5,6,7,8,5},
                      {12,9,10,11,12,9},
                      {16,13,14,15,16,13},
                      {0,1,2,3,4,0}};

/*
if(rank==0){
//printer

double** out_grid= (double**)allocMat(M,sizeof(double));
int recive_size = (size-2)*(size-2);
double* temp = (double*)malloc(recive_size*sizeof(double));
    for(double timestep=0; timestep<T; timestep+=dt)
    {
        while(1){
            MPI_Recv(temp, recive_size, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,&Stat);
            //use Stat.MPI_RANK to locate in array
            

            
        }
    MPI_Barrier(MPI_COMM_WORLD); //wait till next timestep
    }
free(temp);
freeMat(out_grid,M);
}

*/
if(rank!=0){
//find where you are in the topology
int my_x_place;
int my_y_place;
for(int i=1;i<5;i++){
    for(int j=1;j<5;j++){
        if (rank_layout[i][j]==rank){
            my_x_place=i;
            my_y_place=j;
        }
    }
}

//allocMat actually calls calloc, thus is zeroed
double** grid_ptr= (double**)allocMat(size+2,sizeof(double));
double** next = (double**)allocMat(size+2,sizeof(double));
double** prior = (double**)allocMat(size+2,sizeof(double));

double** temp;
    for(double timestep=0; timestep<T; timestep+=dt)
    {
            for(int i=1;i<size+1;i++)
            {
                for(int j=1;j<size+1;j++)
                {
                    //first if statment defines our osilator
                    //TODO: find a better what to define osilators

                    next[i][j]=discrete_wave_eq(grid_ptr,prior,dx,dy,dt,i,j);
                }
            }

            //exchange with other ranks
            //up
            MPI_Sendrecv(next[1],size,MPI_DOUBLE,rank_layout[my_x_place-1][my_y_place],0,
                        next[0],size,MPI_DOUBLE,rank_layout[my_x_place+1][my_y_place],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            //down
            MPI_Sendrecv(next[size-2],size,MPI_DOUBLE,rank_layout[my_x_place+1][my_y_place],0,
                        next[size-1],size,MPI_DOUBLE,rank_layout[my_x_place-1][my_y_place],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            //flip mat
            flipMat_d(next, size+2);
            //left
            MPI_Sendrecv(next[1],size,MPI_DOUBLE,rank_layout[my_x_place-1][my_y_place],0,
                        next[0],size,MPI_DOUBLE,rank_layout[my_x_place+1][my_y_place],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            //right
            MPI_Sendrecv(next[size-2],size,MPI_DOUBLE,rank_layout[my_x_place+1][my_y_place],0,
                        next[size-1],size,MPI_DOUBLE,rank_layout[my_x_place-1][my_y_place],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            //flip back
            flipMat_d(next, size+2);

            //send_to_printer(grid_ptr, size+2, MPI_COMM_WORLD);
            clear(prior,size+2);
            //ptr swap for next time step
            temp = next;
            next=prior;
            prior = grid_ptr;
            grid_ptr = temp;
            //MPI_Barrier(MPI_COMM_WORLD);//block till next timestep
    }
//cleanup
freeMat(grid_ptr,size+2);
freeMat(next,size+2);
freeMat(prior,size+2);
}

MPI_Finalize();
return 0;
}


double discrete_wave_eq(double** u , double** um, double dx, double dy, double dt, int i, int j){
    //bounds checking is done in loop, this will go out of bounds if you let it
    //printf("here\n");
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


void send_to_printer(double**grid, int size, MPI_Comm comm){

double* sendbuf = (double*)malloc((size-2)*(size-2)*sizeof(double));
int t=0;
            for(int i=1;i<size;i++)
            {
                for(int j=1;j<size;j++)
                {
                sendbuf[t] = grid[i][j];
                t++;
            }}

int sendsz = (size-2)*(size-2);
MPI_Send(sendbuf,sendsz, MPI_DOUBLE, 0, 0, comm);
free(sendbuf);
}



