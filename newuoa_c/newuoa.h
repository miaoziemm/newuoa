#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>

// matrix
float *matrix1f(int m);
float **matrix2f(int m,int n);
float ***matrix3f(int x,int y,int z);
double *matrix1d(int m);
double **matrix2d(int m,int n);
double ***matrix3d(int x,int y,int z);
double *slice1d(double *W, int start, int end);
double min2d(double x,double y);
double max2d(double x,double y);
double max3d(double x, double y, double z);
double **matrix_reshape2d(int n, double *w);
//calfun
double calfun(int N, double *X, double F, double **Y);

//trsapp
int trsapp(int N, int NPT, double *XOPT, double **XPT, double *GQ, double *HQ, double *PQ, double DELTA, double *STEP, double *D, double *G, double *HD, double *HS, double CRVMIN);

//biglag
double biglag(int N, int NPT, double *XOPT, double **XPT, double **BMAT, double **ZMAT, int IDZ, int NDIM, int KNEW, double DELTA, double *D, double ALPHA, double *HCOL, double *GC, double *GD, double *S, double *W);

//bigden
int bigden(int N, int NPT, double *XOPT, double **XPT, double **BMAT, double **ZMAT, int IDZ, int NDIM, int KOPT, int KNEW, double *D, double *W, double *VLAG, double BETA, double *S, int WWVEC, int WPROD);

//update
int update(int N, int NPT, double **BMAT, double **ZMAT, int IDZ, int NDIM, double *VLAG, double BETA, int KNEW, double *W);

//newuob
int newuob(double *WW, int N, int NPT, double *X, double RHOBEG, double RHOEND, int IPRINT, int MAXFUN, double *XBASE, double *XOPT, double *XNEW, double *WXPT, double *FVAL, double *GQ, double *HQ, double *PQ, double *WBMAT, double *WZMAT, int NDIM, double *D, double *VLAG, double *W);

//newuoa
int newuoa(int N, int NPT,double *X, double RHOBEG, double RHOEND, int IPRINT, int MAXFUN, double *W);

