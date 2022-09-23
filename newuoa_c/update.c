#include"newuoa.h"

int update(int N, int NPT, double **BMAT, double **ZMAT, int IDZ, int NDIM, double *VLAG, double BETA, int KNEW, double *W){
    /*The arrays BMAT and ZMAT with IDZ are updated, in order to shift the
     interpolation point that has index KNEW. On entry, VLAG contains the
     components of the vector Theta*Wcheck+e_b of the updating formula
     (6.11), and BETA holds the value of the parameter that has this name.
     The vector W is used for working space.*/

    // Set some constants.
    double ONE=1.0, ZERO=0.0;
    int NPTM=NPT-N-1;
    int i,j,k;

    // Apply the rotations that put zeros in the KNEW-th row of ZMAT.
    int JL=1;
    double TEMP,TEMPA,TEMPB;
    for(j=1;j<NPTM;j++){
        if(j==IDZ){JL=IDZ+1;}
        else if(ZMAT[KNEW][j]!=ZERO){
            TEMP=sqrt(ZMAT[KNEW][JL-1]*ZMAT[KNEW][JL-1]+ZMAT[KNEW][j]*ZMAT[KNEW][j]);
            TEMPA=ZMAT[KNEW][JL-1]/TEMP;
            TEMPB=ZMAT[KNEW][j]/TEMP;
            for(i=0;i<NPT;i++){
                TEMP=TEMPA*ZMAT[i][JL-1]+TEMPB*ZMAT[i][j];
                ZMAT[i][j]=TEMPA*ZMAT[i][j]-TEMPB*ZMAT[i][JL-1];
                ZMAT[i][JL-1]=TEMP;
            }
            ZMAT[KNEW][j]=ZERO;
        }
        continue;
    }

    // Put the first NPT components of the KNEW-th column of HLAG into W,
    // and calculate the parameters of the updating formula.
    TEMPA=ZMAT[KNEW][0];
    if(IDZ>=1){TEMPA=-TEMPA;}
    if(JL>1){TEMPB=ZMAT[KNEW][JL-1];}
    for(i=0;i<NPT;i++){
        W[i]=TEMPA*ZMAT[i][0];
        if(JL>1){W[i]=W[i]+TEMPB*ZMAT[i][JL-1];}
        continue;
    }
    
    double ALPHA=W[KNEW];
    double TAU=VLAG[KNEW];
    double TAUSQ=TAU*TAU;
    double DENOM=ALPHA*BETA+TAUSQ;
    VLAG[KNEW]=VLAG[KNEW]-ONE;

    // Complete the updating of ZMAT when there is only one nonzero element
    // in the KNEW-th row of the new matrix ZMAT, but, if IFLAG is set to one,
    // then the first column of ZMAT will be exchanged with another one later.
    int IFLAG=0, JA, JB;
    double SCALA,SCALB;
    if(JL==1){
        TEMP=sqrt(fabs(DENOM));
        TEMPB=TEMPA/TEMP;
        TEMPA=TAU/TEMP;
        for(i=0;i<NPT;i++){ZMAT[i][0]=TEMPA*ZMAT[i][0]-TEMPB*VLAG[i];}
        if(IDZ==0 && TEMP<ZERO){IDZ=1;}
        if(IDZ>=1 && TEMP>=ZERO){IFLAG=1;} 
    }
    else{
        // Complete the updating of ZMAT in the alternative case.
        JA=1;
        if(BETA>=ZERO){JA=JL;}
        JB=JL+1-JA;
        TEMP=ZMAT[KNEW][JB-1]/DENOM;
        TEMPA=TEMP*BETA;
        TEMPB=TEMP*TAU;
        TEMP=ZMAT[KNEW][JA-1];
        SCALA=ONE/sqrt(fabs(BETA)*TEMP*TEMP+TAUSQ);
        SCALB=SCALA*sqrt(fabs(DENOM));
        for(i=0;i<NPT;i++){
            ZMAT[i][JA-1]=SCALA*(TAU*ZMAT[i][JA-1]-TEMP*VLAG[i]);
            ZMAT[i][JB-1]=SCALB*(ZMAT[i][JB-1]-TEMPA*W[i]-TEMPB*VLAG[i]);
        }
        if(DENOM<=ZERO){
            if(BETA<ZERO){IDZ=IDZ+1;}
            if(BETA>=ZERO){IFLAG=1;}
        }
       
    }

    // IDZ is reduced in the following case, and usually the first column
    // of ZMAT is exchanged with a later one.
    if(IFLAG==1){
        IDZ=IDZ-1;
        for(i=0;i<NPT;i++){
            TEMP=ZMAT[i][0];
            ZMAT[i][0]=ZMAT[i][IDZ];
            ZMAT[i][IDZ]=TEMP;
        }
    }

    // Finally, update the matrix BMAT.
    int JP;
    for(j=0;j<N;j++){
        JP=NPT+j;
        W[JP]=BMAT[KNEW][j];
        TEMPA=(ALPHA*VLAG[JP]-TAU*W[JP])/DENOM;
        TEMPB=(-BETA*W[JP]-TAU*VLAG[JP])/DENOM;
        for(i=0;i<=JP;i++){
            BMAT[i][j]=BMAT[i][j]+TEMPA*VLAG[i]+TEMPB*W[i];
            if(i+1>NPT){BMAT[JP][i-NPT]=BMAT[i][j];}
        }
    }
    return 0;
}
