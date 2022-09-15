#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// matrix
float *matrix1f(int m);
float **matrix2f(int m,int n);
float ***matrix3f(int x,int y,int z);
double *matrix1d(int m);
double **matrix2d(int m,int n);
double ***matrix3d(int x,int y,int z);

//trsapp
int trsapp(int N, int NPT, double *XOPT, double **XPT, double *GQ, double *HQ, double *PQ, double DELTA, double *STEP, double *D, double *G, double *HD, double *HS, double CRVMIN);