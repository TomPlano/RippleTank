#include "Utils.h"
#include <stdlib.h> 
#include <stdio.h>
#include <omp.h>
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
    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
        #pragma omp for schedule(dynamic)
        for (int i =0; i< n; i++)
        {
            resvec[i]= dot_prod(vec,mat[i], n);
        } 
    }
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
    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
        #pragma omp for schedule(dynamic)
        for (int i=0;i<size; i++)
        {
            for (int j=0;j<size; j++)
            {
                resMat[i][j] = dot_prod(a[i],b[j],size);
            }
        }
    }
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
    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
        #pragma omp for schedule(dynamic)
        for (int i =0; i< n; i++)
        {
            resvec[i]= dot_prod(vec,mat[i], n);
        } 
    }
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
    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
        #pragma omp for schedule(dynamic)
        for (int i=0;i<size; i++)
        {
            for (int j=0;j<size; j++)
            {
                resMat[i][j] = dot_prod(a[i],b[j],size);
            }
        }
    }
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
    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
        #pragma omp for schedule(dynamic)
        for (int i =0; i< n; i++)
        {
            resvec[i]= dot_prod(vec,mat[i], n);
        } 
    }
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
    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
        #pragma omp for schedule(dynamic)
        for (int i=0;i<size; i++)
        {
            for (int j=0;j<size; j++)
            {
                resMat[i][j] = dot_prod(a[i],b[j],size);
            }
        }
    }
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





