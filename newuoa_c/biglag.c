#include<newuoa.h>

int biglag(int N, int NPT, double *XOPT, double **XPT, double **BMAT, double **ZMAT, int IDZ, int NDIM, int KNEW, double DELTA, double *D, double ALPHA, double *HCOL, double *GC, double *GD, double *S, double *W){
    
    /* N is the number of variables.
     NPT is the number of interpolation equations.
     XOPT is the best interpolation point so far.
     XPT contains the coordinates of the current interpolation points.
     BMAT provides the last N columns of H.
     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
     NDIM is the first dimension of BMAT and has the value NPT+N.
     KNEW is the index of the interpolation point that is going to be moved.
     DELTA is the current trust region bound.
     D will be set to the step from XOPT to the new point.
     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
     HCOL, GC, GD, S and W will be used for working space.

     The step D is calculated in a way that attempts to maximize the modulus
     of LFUNC(XOPT+D), subject to the bound ||D|| .LE. DELTA, where LFUNC is
     the KNEW-th Lagrange function.

     Set some constants. */
    
    int i,j,k;
    double HALF=0.5,ONE=1.0,ZERO=0.0;
    double TWOPI=8.0*atan(ONE);
    double DELSQ=DELTA*DELTA;
    int NPTM=NPT-N-1;

    // Set the first NPT components of HCOL to the leading elements of the KNEW-th column of H.
    int ITREC=0;
    double TEMP;
    for(i=0;i<NPT;i++){
        HCOL[i]=ZERO;
    }
    for(j=0;j<NPTM;j++){
        TEMP=ZMAT[KNEW][j];
        if(j<IDZ){TEMP=-TEMP;}
        for(k=0;k<NPT;k++){
            HCOL[k]=HCOL[k]+TEMP*ZMAT[k][j];
        }
    }
    ALPHA=HCOL[KNEW];

    // Set the unscaled initial direction D. Form the gradient of LFUNC at
    // XOPT, and multiply D by the second derivative matrix of LFUNC.
    double DD=ZERO, SUM;
    for(i=0;i<N;i++){
        D[i]=XPT[KNEW][i]-XOPT[i];
        GC[i]=BMAT[KNEW][i];
        GD[i]=ZERO;
        DD=DD+D[i]*D[i];
    } 
    for(k=0;k<NPT;k++){
        TEMP=ZERO;
        SUM=ZERO;
        for(j=0;j<N;j++){
            TEMP=TEMP+XPT[k][j]*XOPT[j];
            SUM=SUM+XPT[k][j]*D[j];
        }
        TEMP=HCOL[k]*TEMP;
        SUM=HCOL[k]*SUM;
        for(i=0;i<N;i++){
            GC[i]=GC[i]+TEMP*XPT[k][i];
            GD[i]=GD[i]+SUM*XPT[k][i];
        }
    }   

    // Scale D and GD, with a sign change if required. Set S to another
    // vector in the initial two dimensional subspace.
    double GG=ZERO,SP=ZERO,DHD=ZERO;
    for(i=0;i<N;i++){
        GG=GG+GC[i]*GC[i];
        SP=SP+D[i]*GC[i];
        DHD=DHD+D[i]*GD[i];
    }
    double SCALE=DELTA/sqrt(DD);
    if(SP*DHD<ZERO){SCALE=-SCALE;}
    TEMP=ZERO;
    if(SP*SP>0.99*DD*GG){TEMP=ONE;}
    double TAU=SCALE*(abs(SP)+HALF*SCALE*abs(DHD));
    if(GG*DELSQ<0.01*TAU*TAU){TEMP=ONE;}
    for(i=0;i<N;i++){
        D[i]=SCALE*D[i];
        GD[i]=SCALE*GD[i];
        S[i]=GC[i]+TEMP*GD[i];
    }

    // Begin the iteration by overwriting S with a vector that has the
    // required length and direction, except that termination occurs if
    // the given D and S are nearly parallel.
h:    
    ITREC=ITREC+1;
    DD=ZERO;
    SP=ZERO;
    double SS=ZERO;
    for(i=0;i<N;i++){
        DD=DD+D[i]*D[i];
        SP=SP+D[i]*S[i];
        SS=SS+S[i]*S[i];
    }
    TEMP=DD*SS-SP*SP;
    if(TEMP<=1.0e-8*DD*SS){goto p;}
    double DENOM=sqrt(TEMP);
    for(i=0;i<N;i++){
        S[i]=(DD*S[i]-SP*D[i])/DENOM;
        W[i]=ZERO;
    }

    // Calculate the coefficients of the objective function on the circle,
    // beginning with the multiplication of S by the second derivative matrix.
    for(k=0;k<NPT;k++){
        SUM=ZERO;
        for(j=0;j<N;j++){
            SUM=SUM+XPT[k][j]*S[j];
        }
        SUM=HCOL[k]*SUM;
        for(i=0;i<N;i++){
            W[i]=W[i]+SUM*XPT[k][i];
        }
    }
    double CF1=ZERO,CF2=ZERO,CF3=ZERO,CF4=ZERO,CF5=ZERO;
    for(i=0;i<N;i++){
        CF1=CF1+S[i]*W[i];
        CF2=CF2+D[i]*GC[i];
        CF3=CF3+S[i]*GC[i];
        CF4=CF4+D[i]*GD[i];
        CF5=CF5+S[i]*GD[i];
    }
    CF1=HALF*CF1;
    CF4=HALF*CF4-CF1;

    // Seek the value of the angle that maximizes the modulus of TAU.
    double TAUBEG=CF1+CF2+CF3;
    double TAUMAX=TAUBEG;
    double TAUOLD=TAUBEG;
    double ANGLE,CTH,STH,TEMPA,TEMPB;
    int ISAVE=0,IU=49;
    TEMP=TWOPI/(double)(IU+1);
    for(i=1;i<=IU;i++){
        ANGLE=(double)i*TEMP;
        CTH=cos(ANGLE);
        STH=sin(ANGLE);
        TAU=CF1+(CF2+CF4*CTH)+(CF3+CF5*CTH)*STH;
        if(abs(TAU)>abs(TAUMAX)){
            TAUMAX=TAU;
            ISAVE=i;
            TEMPA=TAUOLD;
        }
        else if(i == ISAVE+1){
            TEMPB=TAU;
        }
        TAUOLD=TAU;
    }
    if(ISAVE==0){TEMPA=TAU;}
    if(ISAVE==IU){TEMPB=TAUBEG;}
    double STEP=ZERO;
    if(TEMPA!=TEMPB){
        TEMPA=TEMPA-TAUMAX;
        TEMPB=TEMPB-TAUMAX;
        STEP=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB);
    }
    ANGLE=TEMP*((double)ISAVE+STEP);

    // Calculate the new D and GD. Then test for convergence.
    CTH=cos(ANGLE);
    STH=sin(ANGLE);
    TAU=CF1+(CF2+CF4*CTH)*CTH+(CF3+CF5*CTH)*STH;
    for(i=0;i<N;i++){
        D[i]=CTH*D[i]+STH*S[i];
        GD[i]=CTH*GD[i]+STH*W[i];
        S[i]=GC[i]+GD[i];
    }
    if(abs(TAU)<=1.1*abs(TAUBEG)){goto p;}
    if(ITREC < N){goto h;}
p:    
    return 0;
}