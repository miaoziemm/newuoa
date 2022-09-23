#include "newuoa.h"

double calfun(int N, double *X, double F, double **Y){

    
    int i,j,k,NP;
    for(j=0;j<N;j++){
        Y[0][j]=1.0;
        Y[1][j]=2.0*X[j]-1.0;
    }
    for(i=2;i<=N;i++){
        for(j=1;j<=N;j++){
            Y[i+1-1][j-1]=2.0*Y[1][j-1]*Y[i-1][j-1]-Y[i-1-1][j-1];
            
        }
    }
    F=0.0;
    NP=N+1;
    int IW=1,t;
    double SUM;
    for(i=0;i<NP;i++){
        SUM=0.0;
        for(j=0;j<N;j++){SUM=SUM+Y[i][j];}
        SUM=SUM/(double)N;
        t=i+1;
        if(IW>0){SUM=SUM+1.0/(double)(t*t-2*t);}
        IW=-IW;
        F=F+SUM*SUM;
    }

    return F;
}