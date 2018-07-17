#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <math.h>
#include "Utils.h"

void pgrid(int** grid, int size);
void clear(int** grid,int size);
int discrete_wave_eq(int** grid_t , int** grid_t_minus_one, int size, int x_pos, int y_pos);



int main(int argc, char** argv)
{

if(argc!=3)
{
    printf("./ripple.x <odd size> <itterations>\n");
    return 0;
}


int size = atoi(argv[1]);
int itter = atoi(argv[2]);

int** grid_ptr= (int**)allocMat(size,sizeof(int));//sets all to zero
int** next = (int**)allocMat(size,sizeof(int));
int** prior = (int**)allocMat(size,sizeof(int));

int** temp;


//init t=0 and t=-1
grid_ptr[4][4]=10000;
prior[4][4]=10000;

#pragma omp parallel
{
    for(int timestep=0; timestep<itter; timestep ++)
    {
            
            #pragma omp for schedule(dynamic) 
            for(int i=1;i<size-1;i++)
            {
                    for(int j=1;j<size-1;j++)
                    {
                        next[i][j]=discrete_wave_eq(grid_ptr,prior,size,i,j);
                    }
            }
        #pragma omp master
        {
            //print out cuttent ittereation
            pgrid(grid_ptr,size);
            printf("\n");
            //
            clear(prior,size);
            temp = next;
            next=prior;
            prior = grid_ptr;
            grid_ptr = temp;
        }
        #pragma omp barrier
    }
}


return 0;
}


int discrete_wave_eq(int** grid_t , int** grid_t_minus_one, int size, int x_pos, int y_pos)
{
return 2*grid_t[x_pos][y_pos]-grid_t_minus_one[x_pos][y_pos]+floor((grid_t[x_pos+1][y_pos]+grid_t[x_pos-1][y_pos]+grid_t[x_pos][y_pos+1]+grid_t[x_pos][y_pos-1]-4*grid_t[x_pos][y_pos])/4);
}



void pgrid(int** grid, int size)
{
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<size;j++)
        {
            printf("%i\t",grid[i][j]);
        }
        printf("\n");
    } 
}


void clear(int** grid,int size)
{
  for (int i = 0; i<size-1; i++)
  memset(grid[i],0,size);
  
}
