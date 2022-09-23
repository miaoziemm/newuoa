#include "newuoa.h"

int bigden(int N, int NPT, double *XOPT, double **XPT, double **BMAT, double **ZMAT, int IDZ, int NDIM, int KOPT, int KNEW, double *D, double *W, double *VLAG, double BETA, double *S, int WWVEC, int WPROD){
    /*N is the number of variables.
     NPT is the number of interpolation equations.
     XOPT is the best interpolation point so far.
     XPT contains the coordinates of the current interpolation points.
     BMAT provides the last N columns of H.
     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
     NDIM is the first dimension of BMAT and has the value NPT+N.
     KOPT is the index of the optimal interpolation point.
     KNEW is the index of the interpolation point that is going to be moved.
     D will be set to the step from XOPT to the new point, and on entry it
       should be the D that was calculated by the last call of BIGLAG. The
       length of the initial D provides a trust region bound on the final D.
     W will be set to Wcheck for the final choice of D.
     VLAG will be set to Theta*Wcheck+e_b for the final choice of D.
     BETA will be set to the value that will occur in the updating formula
       when the KNEW-th interpolation point is moved to its new position.
     S, WVEC, PROD and the private arrays DEN, DENEX and PAR will be used
       for working space.

     D is calculated in a way that should provide a denominator with a large
     modulus in the updating formula when the KNEW-th interpolation point is
     shifted to the new position XOPT+D.*/


    // Set some constants.
    
    double PROD[NDIM][WPROD];
    double WVEC[NDIM][WWVEC];

    double *DEN,*DENEX,*PAR;
    DEN=matrix1d(9);DENEX=matrix1d(9);PAR=matrix1d(9);
    int i,j,k,jc;
    double HALF=0.5,ONE=1.0,QUART=0.25,TWO=2.0,ZERO=0.0,TWOPI=8.0*atan(ONE);
    int NPTM=NPT-N-1;

    // Store the first NPT elements of the KNEW-th column of H in W(N+1) to W(N+NPT).
    double TEMP;
    for(k=0;k<NPT;k++){W[N+k]=ZERO;}
    for(j=0;j<NPTM;j++){
        TEMP=ZMAT[KNEW][j];
        if(j<IDZ){TEMP=-TEMP;}
        for(k=0;k<NPT;k++){
            W[N+k]=W[N+k]+TEMP*ZMAT[k][j];
        }
    }
    double ALPHA=W[N+KNEW];

    // The initial search direction D is taken from the last call of BIGLAG,
    // and the initial S is set below, usually to the direction from X_OPT
    // to X_KNEW, but a different direction to an interpolation point may
    // be chosen, in order to prevent S from being nearly parallel to D.
    double DD=ZERO,DS=ZERO,SS=ZERO;
    double XOPTSQ=ZERO;
    for(i=0;i<N;i++){
        DD=DD+D[i]*D[i];
        S[i]=XPT[KNEW][i]-XOPT[i];
        DS=DS+D[i]*S[i];
        SS=SS+S[i]*S[i];
        XOPTSQ=XOPTSQ+XOPT[i]*XOPT[i];
    }
    int KSAV;
    double DSTEMP,SSTEMP,DTEST,DIFF;
    if(DS*DS>0.99*DD*SS){
        KSAV=KNEW;
        DTEST=DS*DS/SS;
        for(k=0;k<NPT;k++){
            if(k!=KOPT){
                DSTEMP=ZERO;
                SSTEMP=ZERO;
                for(i=0;i<N;i++){
                    DIFF=XPT[k][i]-XOPT[i];
                    DSTEMP=DSTEMP+D[i]*DIFF;
                    SSTEMP=SSTEMP+DIFF*DIFF;
                }
                if(DSTEMP*DSTEMP/SSTEMP<DTEST){
                    KSAV=k;
                    DTEST=DSTEMP*DSTEMP/SSTEMP;
                    DS=DSTEMP;
                    SS=SSTEMP;
                }
            }
            continue;
        }
        for(i=0;i<N;i++){S[i]=XPT[KSAV][i]-XOPT[i];}
    }
    double SSDEN=DD*SS-DS*DS,DENSAV=ZERO;
    int ITERC=0;

    // Begin the iteration by overwriting S with a vector that has the
    // required length and direction.
n70:
    ITERC=ITERC+1;
    TEMP=ONE/sqrt(SSDEN);
    double XOPTD=ZERO,XOPTS=ZERO;
    for(i=0;i<N;i++){
        S[i]=TEMP*(DD*S[i]-DS*D[i]);
        XOPTD=XOPTD+XOPT[i]*D[i];
        XOPTS=XOPTS+XOPT[i]*S[i];
    }

    // Set the coefficients of the first two terms of BETA.
    double TEMPA,TEMPB,TEMPC;
    TEMPA=HALF*XOPTD*XOPTD;
    TEMPB=HALF*XOPTS*XOPTS;
    DEN[0]=DD*(XOPTSQ+HALF*DD)+TEMPA+TEMPB;
    DEN[1]=TWO*XOPTD*DD;
    DEN[2]=TWO*XOPTS*DD;
    DEN[3]=TEMPA-TEMPB;
    DEN[4]=XOPTD*XOPTS;
    for(i=5;i<9;i++){DEN[i]=ZERO;}

    // Put the coefficients of Wcheck in WVEC.
    int IP;
    for(k=0;k<NPT;k++){
        TEMPA=ZERO;TEMPB=ZERO;TEMPC=ZERO;
        for(i=0;i<N;i++){
            TEMPA=TEMPA+XPT[k][i]*D[i];
            TEMPB=TEMPB+XPT[k][i]*S[i];
            TEMPC=TEMPC+XPT[k][i]*XOPT[i];
        }
        WVEC[k][0]=QUART*(TEMPA*TEMPA+TEMPB*TEMPB);
        WVEC[k][1]=TEMPA*TEMPC;
        WVEC[k][2]=TEMPB*TEMPC;
        WVEC[k][3]=QUART*(TEMPA*TEMPA-TEMPB*TEMPB);
        WVEC[k][4]=HALF*TEMPA*TEMPB;    
    }
    
    for(i=0;i<N;i++){
        IP=i+NPT;
        WVEC[IP][0]=ZERO;
        WVEC[IP][1]=D[i];
        WVEC[IP][2]=S[i];
        WVEC[IP][3]=ZERO;
        WVEC[IP][4]=ZERO;
    }
    

    // Put the coefficents of THETA*Wcheck in PROD.
    int NW;
    double SUM;
   
    for(jc=0;jc<5;jc++){
        NW=NPT;
        if(jc==1||jc==2){NW=NDIM;}
        for(k=0;k<NPT;k++){PROD[k][jc]=ZERO;}
         
        for(j=0;j<=NPTM-1;j++){
            SUM=ZERO;
            for(k=0;k<NPT;k++){SUM=SUM+ZMAT[k][j]*WVEC[k][jc];}
            if(j<=IDZ-1){SUM=-SUM;}
            for(k=0;k<NPT;k++){PROD[k][jc]=PROD[k][jc]+SUM*ZMAT[k][j];}
        }
        if(NW==NDIM){
            for(k=0;k<NPT;k++){
                SUM=ZERO;
                for(j=0;j<N;j++){
                    SUM=SUM+BMAT[k][j]*WVEC[NPT+j][jc];
                }
                PROD[k][jc]=PROD[k][jc]+SUM;
            }
        }

        for(j=0;j<N;j++){
            SUM=ZERO;
            for(i=0;i<NW;i++){
                SUM=SUM+BMAT[i][j]*WVEC[i][jc];
            }
            PROD[NPT+j][jc]=SUM;
        }
    
    }

    // Include in DEN the part of BETA that depends on THETA.
    for(k=0;k<NDIM;k++){
        SUM=ZERO;
        for(i=0;i<5;i++){
            PAR[i]=HALF*PROD[k][i]*WVEC[k][i];
            SUM=SUM+PAR[i];
        }
        DEN[0]=DEN[0]-PAR[0]-SUM;
        TEMPA=PROD[k][0]*WVEC[k][1]+PROD[k][1]*WVEC[k][0];
        TEMPB=PROD[k][1]*WVEC[k][3]+PROD[k][3]*WVEC[k][1];
        TEMPC=PROD[k][2]*WVEC[k][4]+PROD[k][4]*WVEC[k][2];
        DEN[1]=DEN[1]-TEMPA-HALF*(TEMPB+TEMPC);
        DEN[5]=DEN[5]-HALF*(TEMPB-TEMPC);
        TEMPA=PROD[k][0]*WVEC[k][2]+PROD[k][2]*WVEC[k][0];
        TEMPB=PROD[k][1]*WVEC[k][4]+PROD[k][4]*WVEC[k][1];
        TEMPC=PROD[k][2]*WVEC[k][3]+PROD[k][3]*WVEC[k][2];
        DEN[2]=DEN[2]-TEMPA-HALF*(TEMPB-TEMPC);
        DEN[6]=DEN[6]-HALF*(TEMPB+TEMPC);
        TEMPA=PROD[k][0]*WVEC[k][3]+PROD[k][3]*WVEC[k][0];
        DEN[3]=DEN[3]-TEMPA-PAR[1]+PAR[2];
        TEMPA=PROD[k][0]*WVEC[k][4]+PROD[k][4]*WVEC[k][0];
        TEMPB=PROD[k][1]*WVEC[k][2]+PROD[k][2]*WVEC[k][1];
        DEN[4]=DEN[4]-TEMPA-HALF*TEMPB;
        DEN[7]=DEN[7]-PAR[3]+PAR[4];
        TEMPA=PROD[k][3]*WVEC[k][4]+PROD[k][4]*WVEC[k][3];
        DEN[8]=DEN[8]-HALF*TEMPA;
    }

    // Extend DEN so that it holds all the coefficients of DENOM.
    SUM=ZERO;
    for(i=0;i<5;i++){
        PAR[i]=HALF*PROD[KNEW][i]*PROD[KNEW][i];
        SUM=SUM+PAR[i];
    }
    DENEX[0]=ALPHA*DEN[0]+PAR[0]+SUM;
    TEMPA=TWO*PROD[KNEW][0]*PROD[KNEW][1];
    TEMPB=PROD[KNEW][1]*PROD[KNEW][3];
    TEMPC=PROD[KNEW][2]*PROD[KNEW][4];
    DENEX[1]=ALPHA*DEN[1]+TEMPA+TEMPB+TEMPC;
    DENEX[5]=ALPHA*DEN[5]+TEMPB-TEMPC;
    TEMPA=TWO*PROD[KNEW][0]*PROD[KNEW][2];
    TEMPB=PROD[KNEW][1]*PROD[KNEW][4];
    TEMPC=PROD[KNEW][2]*PROD[KNEW][3];
    DENEX[2]=ALPHA*DEN[2]+TEMPA+TEMPB-TEMPC;
    DENEX[6]=ALPHA*DEN[6]+TEMPB+TEMPC;
    TEMPA=TWO*PROD[KNEW][0]*PROD[KNEW][3];
    DENEX[3]=ALPHA*DEN[3]+TEMPA+PAR[1]-PAR[2];
    TEMPA=TWO*PROD[KNEW][0]*PROD[KNEW][4];
    DENEX[4]=ALPHA*DEN[4]+TEMPA+PROD[KNEW][1]*PROD[KNEW][2];
    DENEX[7]=ALPHA*DEN[7]+PAR[3]-PAR[4];
    DENEX[8]=ALPHA*DEN[8]+PROD[KNEW][3]*PROD[KNEW][4];

    // Seek the value of the angle that maximizes the modulus of DENOM.
    SUM=DENEX[0]+DENEX[1]+DENEX[3]+DENEX[5]+DENEX[7];
    double DENOLD=SUM,DENMAX=SUM,ANGLE,SUMOLD;
    int ISAVE=0,IU=49;
    TEMP=TWOPI/(double)(IU+1);
    PAR[0]=ONE;
    for(i=1;i<=IU;i++){
        ANGLE=(double)i*TEMP;
        PAR[1]=cos(ANGLE);
        PAR[2]=sin(ANGLE);
        for(j=3;j<8;j=j+2){
            PAR[j]=PAR[1]*PAR[j-2]-PAR[2]*PAR[j-1];
            PAR[j+1]=PAR[1]*PAR[j-1]-PAR[2]*PAR[j-2];
        }
        SUMOLD=SUM;
        SUM=ZERO;
        for(j=0;j<9;j++){SUM=SUM+DENEX[j]*PAR[j];}
        if(fabs(SUM)>fabs(DENMAX)){
            DENMAX=SUM;
            ISAVE=i;
            TEMPA=SUMOLD;
        }
        else if(i==ISAVE+1){TEMPB=SUM;}
        continue;
    }
    if(ISAVE==0){TEMPA=SUM;}
    if(ISAVE==IU){TEMPB=DENOLD;}
    double STEP=ZERO;
    if(TEMPA!=TEMPB){
        TEMPA=TEMPA-DENMAX;
        TEMPB=TEMPB-DENMAX;
        STEP=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB);
    }
    ANGLE=TEMP*((double)ISAVE+STEP);

    // Calculate the new parameters of the denominator, the new VLAG vector
    // and the new D. Then test for convergence.
    PAR[1]=cos(ANGLE);
    PAR[2]=sin(ANGLE);
    for(j=3;j<8;j=j+2){
        PAR[j]=PAR[1]*PAR[j-2]-PAR[2]*PAR[j-1];
        PAR[j+1]=PAR[1]*PAR[j-1]-PAR[2]*PAR[j-2];
    }
    BETA=ZERO;
    DENMAX=ZERO;
    for(j=0;j<9;j++){
        BETA=BETA+DEN[j]*PAR[j];
        DENMAX=DENMAX+DENEX[j]*PAR[j];
    }
    for(k=0;k<NDIM;k++){
        VLAG[k]=ZERO;
        for(j=0;j<5;j++){VLAG[k]=VLAG[k]+PROD[k][j]*PAR[j];}
    }
    double TAU=VLAG[KNEW];
    DD=ZERO;
    TEMPA=ZERO;
    TEMPB=ZERO;
    for(i=0;i<N;i++){
        D[i]=PAR[1]*D[i]+PAR[2]*S[i];
        W[i]=XOPT[i]+D[i];
        DD=DD+D[i]*D[i];
        TEMPA=TEMPA+D[i]*W[i];
        TEMPB=TEMPB+W[i]*W[i];
    }
    if(ITERC>=N){goto n340;}
    if(ITERC>1){DENSAV=DENSAV>DENOLD?DENSAV:DENOLD;}
    if(fabs(DENMAX)<=1.1*fabs(DENSAV)){goto n340;}
    DENSAV=DENMAX;

    // Set S to half the gradient of the denominator with respect to D.
    // Then branch for the next iteration.
    for(i=0;i<N;i++){
        TEMP=TEMPA*XOPT[i]+TEMPB*D[i]-VLAG[NPT+i];
        S[i]=TAU*BMAT[KNEW][i]+ALPHA*TEMP;
    }
    for(k=0;k<NPT;k++){
        SUM=ZERO;
        for(j=0;j<N;j++){SUM=SUM+XPT[k][j]*W[j];}
        TEMP=(TAU*W[N+k]-ALPHA*VLAG[k])*SUM;
        for(i=0;i<N;i++){S[i]=S[i]+TEMP*XPT[k][i];}
    }
    SS=ZERO;
    DS=ZERO;
    for(i=0;i<N;i++){
        SS=SS+S[i]*S[i];
        DS=DS+D[i]*S[i];
    }
    SSDEN=DD*SS-DS*DS;
    if(SSDEN>=1.0e-8*DD*SS){goto n70;}

    // Set the vector W before the RETURN from the subroutine.
n340:
    for(k=0;k<NDIM;k++){
        W[k]=ZERO;
        for(j=0;j<5;j++){W[k]=W[k]+WVEC[k][j]*PAR[j];}
    }
    VLAG[KOPT]=VLAG[KOPT]+ONE;
    
    return 0;
}