#include "newuoa.h"

float *matrix1f(int m)
{
	float *p;
	p=(float *)malloc(m*sizeof(float));
	return p;
}

float **matrix2f(int m,int n)
{
	float **pp;
	int i,j;
	pp=(float **)malloc(m*sizeof(float *));
	for(i=0;i<m;i++)
	{
		pp[i]=(float *)malloc(n*sizeof(float));
	}
return pp;
}//2-D matrix generation

float ***matrix3f(int x,int y,int z)
{
	float ***ppp;
	int i,j;
	ppp=(float ***)malloc(x*sizeof(float **));
	for(i=0;i<x;i++)
	{
		ppp[i]=(float **)malloc(y*sizeof(float *));
	}
	for(i=0;i<x;i++)
	{
		for(j=0;j<y;j++)
		{
			ppp[i][j]=(float *)malloc(z*sizeof(float));
		}
	}
return ppp;
}//3-D matrix generation

double *matrix1d(int m)
{
	double *p;
	p=(double *)malloc(m*sizeof(double));
	return p;
}

double **matrix2d(int m,int n)
{
	double **pp;
	int i,j;
	pp=(double **)malloc(m*sizeof(double *));
	for(i=0;i<m;i++)
	{
		pp[i]=(double *)malloc(n*sizeof(double));
	}
return pp;
}//2-D matrix generation

double ***matrix3d(int x,int y,int z)
{
	double ***ppp;
	int i,j;
	ppp=(double ***)malloc(x*sizeof(double **));
	for(i=0;i<x;i++)
	{
		ppp[i]=(double **)malloc(y*sizeof(double *));
	}
	for(i=0;i<x;i++)
	{
		for(j=0;j<y;j++)
		{
			ppp[i][j]=(double *)malloc(z*sizeof(double));
		}
	}
return ppp;
}//3-D matrix generation