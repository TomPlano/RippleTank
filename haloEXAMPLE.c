#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {

 int numranks, rank;
 MPI_Init(&argc,&argv);
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 MPI_Comm_size(MPI_COMM_WORLD, &numranks);



 int size=5;
 int data[size+1][size];

//itin
if(rank==0){
    for(int i=0;i<=size;i++){
        for(int j=0; j<size;j++){
            data[i][j]=rank;
            if(i==size)data[i][j]=-1;
        }
    }
}

if(rank==1){

    for(int i=0;i<=size;i++){
        for(int j=0; j<size;j++){
            data[i][j]=rank;
            if(i==0)data[i][j]=-1;
        }
    }
}


//halo

if(rank==0){
 MPI_Sendrecv(data[size-1],size,MPI_INT,
                1,0,data[size],size,MPI_INT,
                1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}
if(rank==1){
 MPI_Sendrecv(data[1],size,MPI_INT,
                0,0,data[0],size,MPI_INT,
                0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}

    

//print
if(rank==0){
    printf("Rank %i\n", rank);
    for(int i=0;i<=size;i++){
        for(int j=0; j<size;j++){
            printf("%i ",data[i][j]);
        }
        printf("\n");
    }
}
MPI_Barrier(MPI_COMM_WORLD);
if(rank==1){
    printf("Rank %i\n", rank);
    for(int i=0;i<=size;i++){
        for(int j=0; j<size;j++){
            printf("%i ",data[i][j]);
        }
        printf("\n");
    }
}


 MPI_Finalize();
}
