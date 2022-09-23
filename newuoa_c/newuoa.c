#include "newuoa.h"

int newuoa(int N, int NPT,double *X, double RHOBEG, double RHOEND, int IPRINT, int MAXFUN, double *W){
     
    int i,j,k;
    int NP,NPTM;
    NP=N+1;
    NPTM=NPT-NP;
    if(NPT<N+2 || NPT>((N+2)*NP)/2){
        printf("Return from NEWUOA because NPT is not in the required interval\n");
        goto n20;
    }
    int NDIM=NPT+N;
    int IXB=1;
    int IXO=IXB+N;
    int IXN=IXO+N;
    int IXP=IXN+N;
    int IFV=IXP+N*NPT;
    int IGQ=IFV+NPT;
    int IHQ=IGQ+N;
    int IPQ=IHQ+(N*NP)/2;
    int IBMAT=IPQ+NPT;
    int IZMAT=IBMAT+NDIM*N;
    int ID=IZMAT+NPT*NPTM;
    int IVL=ID+N;
    int IW=IVL+NDIM;

    double *W_IXB=slice1d(W,IXB-1,10001);
    double *W_IXO=slice1d(W,IXO-1,10001);
    double *W_IXN=slice1d(W,IXN-1,10001);
    double *W_IXP=slice1d(W,IXP-1,10001);
    double *W_IFV=slice1d(W,IFV-1,10001);
    double *W_IGQ=slice1d(W,IGQ-1,10001);
    double *W_IHQ=slice1d(W,IHQ-1,10001);
    double *W_IPQ=slice1d(W,IPQ-1,10001);
    double *W_IBMAT=slice1d(W,IBMAT-1,10001);
    double *W_IZMAT=slice1d(W,IZMAT-1,10001);
    double *W_ID=slice1d(W,ID-1,10001);
    double *W_IVL=slice1d(W,IVL-1,10001);
    double *W_IW=slice1d(W,IW-1,10001);

    newuob(W,N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W_IXB,W_IXO,W_IXN,W_IXP,W_IFV,W_IGQ,W_IHQ,W_IPQ,W_IBMAT,W_IZMAT,NDIM,W_ID,W_IVL,W_IW);
n20:
    return 0;
}