//header file for useful stuff for parallel class

#ifndef _TOMUTILS_
#define _TOMSUTILS_


#define assignMat(X,Y) _Generic((X), \
                double**: assignMat_d, \
                float**:  assignMat_f, \
                long double**:  assignMat_ld \
)(X,Y)
            
#define matVec(W,X,Y,Z) _Generic((W), \
                double**: matVec_d, \
                float**:  matVec_f, \
                long double**:  matVec_ld \
)(W,X,Y,Z)

#define dot_prod(X,Y,Z) _Generic((X), \
                double*: dot_prod_d, \
                float*: dot_prod_f, \
                long double*: dot_prod_ld \
)(X,Y,Z)
#define freeMat(X,Y) _Generic((X), \
                double**: freeMat_d, \
                float**: freeMat_f, \
                long double**: freeMat_ld \
)(X,Y)

#define flipMat(X,Y) _Generic((X), \
                double**: flipMat_d, \
                float**: flipMat_f, \
                long double**: flipMat_ld \
)(X,Y)

#define matMult(W,X,Y,Z) _Generic((W), \
                double**: matMult_d, \
                float**:  matMult_f, \
                long double**: matMult_ld \
)(W,X,Y,Z)

#define printMat(X,Y) _Generic((X), \
                double**: printMat_d, \
                float**:  printMat_f, \
                long double**: printMat_ld \
)(X,Y)




void* allocVec(int size, int type_size);
void** allocMat(int size, int type_size);

//functions for double
void assignMat_d(double** mat,int size);
void matVec_d(double** mat, double* vec, double* resvec,int n);
double dot_prod_d(double* a, double* b, int size);
void freeMat_d(double** mat, int size);
void flipMat_d(double** a, int size);
void matMult_d(double** a, double** b, double** resMat, int size);
void printMat_d(double** a, int size);


//functions for floats
void printMat_f(float** a, int size);
void assignMat_f(float** mat,int size);
void matVec_f(float** mat, float* vec, float* resvec,int n);
float dot_prod_f(float* a, float* b, int size);
void freeMat_f(float** mat, int size);
void flipMat_f(float** a, int size);
void matMult_f(float** a, float** b, float** resMat, int size);

//functions for long double
void printMat_ld(long double** a, int size);
void assignMat_ld(long double** mat,int size);
void matVec_ld(long double** mat, long double* vec, long double* resvec,int n);
long double dot_prod_ld(long double* a, long double* b, int size);
void freeMat_ld(long double** mat, int size);
void flipMat_ld(long double** a, int size);
void matMult_ld(long double** a, long double** b, long double** resMat, int size);



#endif
