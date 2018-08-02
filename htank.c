#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include "Utils.h"

void pgrid(double** grid, int size);
void clear(double** grid,int size);
double discrete_wave_eq(double** u, double** um, double dx, double dy, double dt, int i, int j);
void send_to_printer(double**grid, double* sendbuf, int size, MPI_Comm comm);

int main(int argc, char** argv)
{
        int numranks, rank;
        MPI_Status Stat;
        MPI_Init(&argc,&argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &numranks);
        double** grid_ptr; 
        double** next; 
        double** prior;
        double* sendbuf;
        double** out_grid;

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


        if(rank==0) {
//printer

                out_grid= (double**)allocMat(M,sizeof(double));
                int recive_size = (size)*(size);
                double* temp = (double*)malloc(recive_size*sizeof(double));
                for(double timestep=0; timestep<T; timestep+=dt)
                {
                        int recived_count=0;
                        while(recived_count<16) {
                                MPI_Recv(temp, recive_size, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,&Stat);
                                //use Stat.MPI_SOURCE to locate in array
                                //printf("Recived from %i\n",Stat.MPI_SOURCE);
                                recived_count++;
                                //put data where it belongs
                                int xstart = (Stat.MPI_SOURCE-1)%4 * size;
                                int xend = xstart + (size);
                                int t = 0;
                 
                                for(int i = xstart; i < xend; i++){
                                    for(int j = xstart; j < xend; j++){
                                        out_grid[i][j] = temp[t];
                                        t++;
                                    }
                                } 
                        }
                        //printf("=======================\n");
                        pgrid(out_grid,M);
                        //print out actual data
                        printf("\n");
                        printf("\n");
                        MPI_Barrier(MPI_COMM_WORLD); //wait till next timestep
                }
                free(temp);
//freeMat(out_grid,M);
        }


        if(rank!=0) {
//find where you are in the topology
                int my_x_place;
                int my_y_place;
                for(int i=1; i<5; i++) {
                        for(int j=1; j<5; j++) {
                                if (rank_layout[i][j]==rank) {
                                        my_x_place=i;
                                        my_y_place=j;
                                }
                        }
                }

//allocMat actually calls calloc, thus is zeroed
                grid_ptr= (double**)allocMat(size+2,sizeof(double));
                next = (double**)allocMat(size+2,sizeof(double));
                prior = (double**)allocMat(size+2,sizeof(double));

                double** temp;
                int sendsz = (size)*(size);
                sendbuf = (double*)malloc(sendsz*sizeof(double));

                for(double timestep=0; timestep<T; timestep+=dt)
                {
                        for(int i=1; i<size+1; i++)
                        {
                                for(int j=1; j<size+1; j++)
                                {
                                        //first if statment defines our osilator
                                        //TODO: find a better what to define osilators
                                        if (rank==10&&i==j&&i==size){
                                            next[i][j] = 5*sin(2*M_PI*3*timestep);
                                        }
                                        else{
                                           next[i][j]=discrete_wave_eq(grid_ptr,prior,dx,dy,dt,i,j);
                                        }
                                        next[i][j]=rank;
                                }
                        }

                        //exchange with other ranks
                        //up
                        MPI_Sendrecv(next[1],size,MPI_DOUBLE,rank_layout[my_x_place-1][my_y_place],0,
                                     next[size-2],size,MPI_DOUBLE,rank_layout[my_x_place+1][my_y_place],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                        //down
                        MPI_Sendrecv(next[size-2],size,MPI_DOUBLE,rank_layout[my_x_place+1][my_y_place],0,
                                     next[0],size,MPI_DOUBLE,rank_layout[my_x_place-1][my_y_place],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                        //flip mat
                        flipMat_d(next, size+2);//will cause trouble
                        //left
                        MPI_Sendrecv(next[1],size,MPI_DOUBLE,rank_layout[my_x_place-1][my_y_place],0,
                                     next[size-2],size,MPI_DOUBLE,rank_layout[my_x_place+1][my_y_place],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                        //right
                        MPI_Sendrecv(next[size-2],size,MPI_DOUBLE,rank_layout[my_x_place+1][my_y_place],0,
                                     next[0],size,MPI_DOUBLE,rank_layout[my_x_place-1][my_y_place],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                        //flip back
                        flipMat_d(next, size+2);
                        //prep for send (serialize)
                        int t=0;
                        for(int i=1; i<size+1; i++)
                        {
                                for(int j=1; j<size+1; j++)
                                {
                                        sendbuf[t] = next[i][j];
                                        t++;
                                }
                        }
                        //send to printer
                        MPI_Send(sendbuf,sendsz, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

                        clear(prior,size+2);
                        //ptr swap for next time step
                        temp = next;
                        next=prior;
                        prior = grid_ptr;
                        grid_ptr = temp;
                        MPI_Barrier(MPI_COMM_WORLD);//block till next timestep
                }
                //cleanup
                /*
                   for(int i=0;i<size+2;i++){
                    free(next[i]);
                    free(grid_ptr[i]);
                    free(prior[i]);

                   }

                   free(next);
                   free(grid_ptr);
                   free(prior);
                 */

                //free(sendbuf);
        }

        MPI_Finalize();
        return 0;
}


double discrete_wave_eq(double** u, double** um, double dx, double dy, double dt, int i, int j){
        //bounds checking is done in loop, this will go out of bounds if you let it
        //printf("here\n");
        return 2*u[i][j]-um[i][j]+((dt*dt)/(2*dx*dx))*(u[i][j-1]-2*u[i][j]+u[i][j+1])+((dt*dt)/(2*dy*dy))*(u[i-1][j]-2*u[i][j]+u[i+1][j]);
}


//TODO: support non square sufaces at some point
void pgrid(double** grid, int size){
        for(int i=0; i<size; i++) {
                for(int j=0; j<size; j++) {
                        printf("%.*f ",10,grid[i][j]);
                }
                printf("\n");
        }
}


//TODO: support non square sufaces at some point
void clear(double** grid,int size){
        for (int i = 0; i<size-1; i++)
                memset(grid[i],0,size);

}
