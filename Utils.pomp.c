
#include "Utils.c.opari.inc"
#include "Utils.h"
#include <stdlib.h> 
#include <stdio.h>
void* allocVec(int size, int type_size)
{
    return (void*)malloc(size*type_size); 
}

void** allocMat(int size, int type_size )
{
    void** temp = (void**)calloc(size,sizeof(void*));
    int i;
    for(i = 0; i<size; i++){
        temp[i] = (void*)calloc(size, type_size);
    }
    return temp;
}




// functions for double
double dot_prod_d(double* a, double* b, int size){
    int i;
    double total = 0;
    for(i=0; i<size; i++)
    {
        total+=a[i]*b[i];
    }
}

void matVec_d(double** mat, double* vec, double* resvec, int n)
{ 
    int num_threads;
    double start = omp_get_wtime(); 
{
  int pomp2_num_threads = omp_get_max_threads();
  int pomp2_if = 1;
  POMP2_Task_handle pomp2_old_task;
  POMP2_Parallel_fork(&pomp2_region_1, pomp2_if, pomp2_num_threads, &pomp2_old_task, pomp2_ctc_1 );
    #pragma omp parallel POMP2_DLIST_00001 firstprivate(pomp2_old_task) if(pomp2_if) num_threads(pomp2_num_threads)
{   POMP2_Parallel_begin( &pomp2_region_1 );
    {
        num_threads = omp_get_num_threads();
{   POMP2_For_enter( &pomp2_region_2, pomp2_ctc_2  );
        #pragma omp for schedule(dynamic) nowait
        for (int i =0; i< n; i++)
        {
            resvec[i]= dot_prod(vec,mat[i], n);
        }
{ POMP2_Task_handle pomp2_old_task;
  POMP2_Implicit_barrier_enter( &pomp2_region_2, &pomp2_old_task );
#pragma omp barrier
  POMP2_Implicit_barrier_exit( &pomp2_region_2, pomp2_old_task ); }
  POMP2_For_exit( &pomp2_region_2 );
 }
          
    }
{ POMP2_Task_handle pomp2_old_task;
  POMP2_Implicit_barrier_enter( &pomp2_region_1, &pomp2_old_task );
#pragma omp barrier
  POMP2_Implicit_barrier_exit( &pomp2_region_1, pomp2_old_task ); }
  POMP2_Parallel_end( &pomp2_region_1 ); }
  POMP2_Parallel_join( &pomp2_region_1, pomp2_old_task ); }
    double end= omp_get_wtime(); 
    printf("%i, ", num_threads);
    printf(" %f\n", end-start);
}

void** allocMat_d(int size, int type_size )
{
    void** temp = (void**)malloc(size*sizeof(void*));
    int i;
    for(i = 0; i<size; i++){
        temp[i] = (void*)calloc(size, type_size);
    }
    return temp;
}

void assignMat_d(double** mat,int size)
{
    int i, j;
    for (int i =0; i< size; i++)
    {
        for(int j=0; j<size; j++)
        {
            if(i==j) mat[i][j] = 2;
            else if (i==j-1||j==i-1) mat[i][j]=1;
        }
    }
}

void freeMat_d(double** mat, int size)
{
    int i;
    for(i = 0; i<size; i++)
    {
        free(mat[i]);
    }
    free(mat);
}

void flipMat_d(double** a, int size)
{   
     for (int i=0;i<size; i++)
    {
        for (int j=i;j<size; j++)
        {
            double temp = a[i][j];
            a[i][j] = a[j][i];
            a[j][i] = temp; 
        }
    }
}

void matMult_d(double** a, double** b, double** resMat, int size)
{
    flipMat(b,size);
    double start = omp_get_wtime(); 
    int num_threads;
{
  int pomp2_num_threads = omp_get_max_threads();
  int pomp2_if = 1;
  POMP2_Task_handle pomp2_old_task;
  POMP2_Parallel_fork(&pomp2_region_3, pomp2_if, pomp2_num_threads, &pomp2_old_task, pomp2_ctc_3 );
    #pragma omp parallel POMP2_DLIST_00003 firstprivate(pomp2_old_task) if(pomp2_if) num_threads(pomp2_num_threads)
{   POMP2_Parallel_begin( &pomp2_region_3 );
    {
        num_threads = omp_get_num_threads();
{   POMP2_For_enter( &pomp2_region_4, pomp2_ctc_4  );
        #pragma omp for schedule(dynamic) nowait
        for (int i=0;i<size; i++)
        {
            for (int j=0;j<size; j++)
            {
                resMat[i][j] = dot_prod(a[i],b[j],size);
            }
        }
{ POMP2_Task_handle pomp2_old_task;
  POMP2_Implicit_barrier_enter( &pomp2_region_4, &pomp2_old_task );
#pragma omp barrier
  POMP2_Implicit_barrier_exit( &pomp2_region_4, pomp2_old_task ); }
  POMP2_For_exit( &pomp2_region_4 );
 }
    }
{ POMP2_Task_handle pomp2_old_task;
  POMP2_Implicit_barrier_enter( &pomp2_region_3, &pomp2_old_task );
#pragma omp barrier
  POMP2_Implicit_barrier_exit( &pomp2_region_3, pomp2_old_task ); }
  POMP2_Parallel_end( &pomp2_region_3 ); }
  POMP2_Parallel_join( &pomp2_region_3, pomp2_old_task ); }
    double end = omp_get_wtime();
    printf("%i, ", num_threads);
    printf(" %f\n", end-start);
    flipMat(b,size);
}

void printMat_d(double** a, int size)
{
    for (int i=0;i<size; i++)
    {
        for (int j=0;j<size; j++)
        {
            printf("%f ", a[i][j]);
        }
        printf("\n");
    }
}



// functions for float
float dot_prod_f(float* a, float* b, int size){
    int i;
    float total = 0;
    for(i=0; i<size; i++)
    {
        total+=a[i]*b[i];
    }
}

void matVec_f(float** mat, float* vec, float* resvec, int n)
{ 
    int num_threads;
    double start = omp_get_wtime(); 
{
  int pomp2_num_threads = omp_get_max_threads();
  int pomp2_if = 1;
  POMP2_Task_handle pomp2_old_task;
  POMP2_Parallel_fork(&pomp2_region_5, pomp2_if, pomp2_num_threads, &pomp2_old_task, pomp2_ctc_5 );
    #pragma omp parallel POMP2_DLIST_00005 firstprivate(pomp2_old_task) if(pomp2_if) num_threads(pomp2_num_threads)
{   POMP2_Parallel_begin( &pomp2_region_5 );
    {
        num_threads = omp_get_num_threads();
{   POMP2_For_enter( &pomp2_region_6, pomp2_ctc_6  );
        #pragma omp for schedule(dynamic) nowait
        for (int i =0; i< n; i++)
        {
            resvec[i]= dot_prod(vec,mat[i], n);
        }
{ POMP2_Task_handle pomp2_old_task;
  POMP2_Implicit_barrier_enter( &pomp2_region_6, &pomp2_old_task );
#pragma omp barrier
  POMP2_Implicit_barrier_exit( &pomp2_region_6, pomp2_old_task ); }
  POMP2_For_exit( &pomp2_region_6 );
 }
          
    }
{ POMP2_Task_handle pomp2_old_task;
  POMP2_Implicit_barrier_enter( &pomp2_region_5, &pomp2_old_task );
#pragma omp barrier
  POMP2_Implicit_barrier_exit( &pomp2_region_5, pomp2_old_task ); }
  POMP2_Parallel_end( &pomp2_region_5 ); }
  POMP2_Parallel_join( &pomp2_region_5, pomp2_old_task ); }
    double end= omp_get_wtime(); 
    printf("%i, ", num_threads);
    printf(" %f\n", end-start);
}

void** allocMat_f(int size, int type_size )
{
    void** temp = (void**)malloc(size*sizeof(void*));
    int i;
    for(i = 0; i<size; i++){
        temp[i] = (void*)calloc(size, type_size);
    }
    return temp;
}

void assignMat_f(float** mat,int size)
{
    int i, j;
    for (int i =0; i< size; i++)
    {
        for(int j=0; j<size; j++)
        {
            if(i==j) mat[i][j] = 2;
            else if (i==j-1||j==i-1) mat[i][j]=1;
        }
    }
}

void freeMat_f(float** mat, int size)
{
    int i;
    for(i = 0; i<size; i++)
    {
        free(mat[i]);
    }
}

void flipMat_f(float** a, int size)
{   
     for (int i=0;i<size; i++)
    {
        for (int j=0;j<size; j++)
        {
            float temp = a[i][j];
            a[i][j] = a[j][i];
            a[j][i] = temp; 
        }
    }
}

void matMult_f(float** a, float** b, float** resMat, int size)
{
    flipMat(b,size);
    double start = omp_get_wtime(); 
    int num_threads;
{
  int pomp2_num_threads = omp_get_max_threads();
  int pomp2_if = 1;
  POMP2_Task_handle pomp2_old_task;
  POMP2_Parallel_fork(&pomp2_region_7, pomp2_if, pomp2_num_threads, &pomp2_old_task, pomp2_ctc_7 );
    #pragma omp parallel POMP2_DLIST_00007 firstprivate(pomp2_old_task) if(pomp2_if) num_threads(pomp2_num_threads)
{   POMP2_Parallel_begin( &pomp2_region_7 );
    {
        num_threads = omp_get_num_threads();
{   POMP2_For_enter( &pomp2_region_8, pomp2_ctc_8  );
        #pragma omp for schedule(dynamic) nowait
        for (int i=0;i<size; i++)
        {
            for (int j=0;j<size; j++)
            {
                resMat[i][j] = dot_prod(a[i],b[j],size);
            }
        }
{ POMP2_Task_handle pomp2_old_task;
  POMP2_Implicit_barrier_enter( &pomp2_region_8, &pomp2_old_task );
#pragma omp barrier
  POMP2_Implicit_barrier_exit( &pomp2_region_8, pomp2_old_task ); }
  POMP2_For_exit( &pomp2_region_8 );
 }
    }
{ POMP2_Task_handle pomp2_old_task;
  POMP2_Implicit_barrier_enter( &pomp2_region_7, &pomp2_old_task );
#pragma omp barrier
  POMP2_Implicit_barrier_exit( &pomp2_region_7, pomp2_old_task ); }
  POMP2_Parallel_end( &pomp2_region_7 ); }
  POMP2_Parallel_join( &pomp2_region_7, pomp2_old_task ); }
    double end = omp_get_wtime();
    printf("%i, ", num_threads);
    printf(" %f\n", end-start);
    flipMat(b,size);
}

void printMat_f(float** a, int size)
{
    for (int i=0;i<size; i++)
    {
        for (int j=0;j<size; j++)
        {
            printf("%f ", a[i][j]);
        }
        printf("\n");
    }
}




// functions for long double
long double dot_prod_ld(long double* a, long double* b, int size){
    int i;
    long double total = 0;
    for(i=0; i<size; i++)
    {
        total+=a[i]*b[i];
    }
}

void matVec_ld(long double** mat, long double* vec, long double* resvec, int n)
{ 
    int num_threads;
    double start = omp_get_wtime(); 
{
  int pomp2_num_threads = omp_get_max_threads();
  int pomp2_if = 1;
  POMP2_Task_handle pomp2_old_task;
  POMP2_Parallel_fork(&pomp2_region_9, pomp2_if, pomp2_num_threads, &pomp2_old_task, pomp2_ctc_9 );
    #pragma omp parallel POMP2_DLIST_00009 firstprivate(pomp2_old_task) if(pomp2_if) num_threads(pomp2_num_threads)
{   POMP2_Parallel_begin( &pomp2_region_9 );
    {
        num_threads = omp_get_num_threads();
{   POMP2_For_enter( &pomp2_region_10, pomp2_ctc_10  );
        #pragma omp for schedule(dynamic) nowait
        for (int i =0; i< n; i++)
        {
            resvec[i]= dot_prod(vec,mat[i], n);
        }
{ POMP2_Task_handle pomp2_old_task;
  POMP2_Implicit_barrier_enter( &pomp2_region_10, &pomp2_old_task );
#pragma omp barrier
  POMP2_Implicit_barrier_exit( &pomp2_region_10, pomp2_old_task ); }
  POMP2_For_exit( &pomp2_region_10 );
 }
          
    }
{ POMP2_Task_handle pomp2_old_task;
  POMP2_Implicit_barrier_enter( &pomp2_region_9, &pomp2_old_task );
#pragma omp barrier
  POMP2_Implicit_barrier_exit( &pomp2_region_9, pomp2_old_task ); }
  POMP2_Parallel_end( &pomp2_region_9 ); }
  POMP2_Parallel_join( &pomp2_region_9, pomp2_old_task ); }
    double end= omp_get_wtime(); 
    printf("%i, ", num_threads);
    printf(" %f\n", end-start);
}

void** allocMat_ld(int size, int type_size )
{
    void** temp = (void**)malloc(size*sizeof(void*));
    int i;
    for(i = 0; i<size; i++){
        temp[i] = (void*)calloc(size, type_size);
    }
    return temp;
}

void assignMat_ld(long double** mat,int size)
{
    int i, j;
    for (int i =0; i< size; i++)
    {
        for(int j=0; j<size; j++)
        {
            if(i==j) mat[i][j] = 2;
            else if (i==j-1||j==i-1) mat[i][j]=1;
        }
    }
}

void freeMat_ld(long double** mat, int size)
{
    int i;
    for(i = 0; i<size; i++)
    {
        free(mat[i]);
    }
}

void flipMat_ld(long double** a, int size)
{   
     for (int i=0;i<size; i++)
    {
        for (int j=0;j<size; j++)
        {
            long double temp = a[i][j];
            a[i][j] = a[j][i];
            a[j][i] = temp; 
        }
    }
}

void matMult_ld(long double** a, long double** b, long double** resMat, int size)
{
    flipMat(b,size);
    double start = omp_get_wtime(); 
    int num_threads;
{
  int pomp2_num_threads = omp_get_max_threads();
  int pomp2_if = 1;
  POMP2_Task_handle pomp2_old_task;
  POMP2_Parallel_fork(&pomp2_region_11, pomp2_if, pomp2_num_threads, &pomp2_old_task, pomp2_ctc_11 );
    #pragma omp parallel POMP2_DLIST_00011 firstprivate(pomp2_old_task) if(pomp2_if) num_threads(pomp2_num_threads)
{   POMP2_Parallel_begin( &pomp2_region_11 );
    {
        num_threads = omp_get_num_threads();
{   POMP2_For_enter( &pomp2_region_12, pomp2_ctc_12  );
        #pragma omp for schedule(dynamic) nowait
        for (int i=0;i<size; i++)
        {
            for (int j=0;j<size; j++)
            {
                resMat[i][j] = dot_prod(a[i],b[j],size);
            }
        }
{ POMP2_Task_handle pomp2_old_task;
  POMP2_Implicit_barrier_enter( &pomp2_region_12, &pomp2_old_task );
#pragma omp barrier
  POMP2_Implicit_barrier_exit( &pomp2_region_12, pomp2_old_task ); }
  POMP2_For_exit( &pomp2_region_12 );
 }
    }
{ POMP2_Task_handle pomp2_old_task;
  POMP2_Implicit_barrier_enter( &pomp2_region_11, &pomp2_old_task );
#pragma omp barrier
  POMP2_Implicit_barrier_exit( &pomp2_region_11, pomp2_old_task ); }
  POMP2_Parallel_end( &pomp2_region_11 ); }
  POMP2_Parallel_join( &pomp2_region_11, pomp2_old_task ); }
    double end = omp_get_wtime();
    printf("%i, ", num_threads);
    printf(" %f\n", end-start);
    flipMat(b,size);
}

void printMat_ld(long double** a, int size)
{
    for (int i=0;i<size; i++)
    {
        for (int j=0;j<size; j++)
        {
            printf("%Lf ", a[i][j]);
        }
        printf("\n");
    }
}





