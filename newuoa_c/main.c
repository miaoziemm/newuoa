#include "newuoa.h"

void main(){
    int i,j,k,N,NPT;
    double *X,*W;
    X=matrix1d(10);
    W=matrix1d(10000);
    int IPRINT=2;
    int MAXFUN=5000;
    double RHOEND=1.0e-6,RHOBEG;
    for(N=2;N<=8;N=N+2){
          
        NPT=2*N+1;
        
        for(i=1;i<=N;i++){X[i-1]=(double)i/(double)(N+1);}
        RHOBEG=0.2*X[0];
        printf("Results with N =%d, and NPT =%d\n",N,NPT);
        
        newuoa(N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W);
        continue;
    }
}