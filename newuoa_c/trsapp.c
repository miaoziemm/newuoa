#include "newuoa.h"

int trsapp(int N, int NPT, double *XOPT, double **XPT, double *GQ, double *HQ, double *PQ, double DELTA, double *STEP, double *D, double *G, double *HD, double *HS, double CRVMIN){
    /*N is the number of variables of a quadratic objective function, Q say.
     The arguments NPT, XOPT, XPT, GQ, HQ and PQ have their usual meanings,
       in order to define the current quadratic model Q.
     DELTA is the trust region radius, and has to be positive.
     STEP will be set to the calculated trial step.
     The arrays D, G, HD and HS will be used for working space.
     CRVMIN will be set to the least curvature of H along the conjugate
       directions that occur, except that it is set to zero if STEP goes
       all the way to the trust region boundary.

     The calculation of STEP begins with the truncated conjugate gradient
     method. If the boundary of the trust region is reached, then further
     changes to STEP may be made, each one being in the 2D space spanned
     by the current STEP and the corresponding gradient of Q. Thus STEP
     should provide a substantial reduction to Q within the trust region. */
    
    
    // Initialization, which includes setting HD to H times XOPT.
    int i,j,k;
    double TEMP;
    double HALF=0.5, ZERO=0.0;
    double TWOPI=8.0*atan(1.0);
    double DELSQ=DELTA*DELTA;
    int ITERC=0, ITERMAX=N, ITERSW=ITERMAX, IH;
    for(i=0;i<N;i++){
        D[i]=XOPT[i];
    }
    goto n170;

    // Prepare for the first line search.
n20:
    double QRED=ZERO;
    double DD=ZERO;
    for(i=0;i<N;i++){
      STEP[i]=ZERO;
      HS[i]=ZERO;
      G[i]=GQ[i]+HD[i];
      D[i]=-G[i];
      DD=DD+pow(D[i],2);
    }
    CRVMIN=ZERO;
    if(DD==ZERO){goto n160;}
    double DS=ZERO;
    double SS=ZERO;
    double GG=DD;
    double GGBEG=GG;

    // Calculate the step to the trust region boundary and the product HD.
n40:
    ITERC=ITERC+1;
    TEMP=DELSQ-SS;
    double BSTEP=TEMP/(DS+sqrt(DS*DS+DD*TEMP));
    goto n170;
n50:
    double DHD=ZERO;
    for(j=0;j<N;j++){DHD=DHD+D[j]*HD[j];}

    // Update CRVMIN and set the step-length ALPHA.
    double ALPHA=BSTEP;
    if(DHD>ZERO){
      TEMP=DHD/DD;
      if(ITERC==1){CRVMIN=TEMP;}
      CRVMIN=CRVMIN<TEMP?CRVMIN:TEMP;
      ALPHA=ALPHA<GG/DHD?ALPHA:GG/DHD;
    }
    double QADD=ALPHA*(GG-HALF*ALPHA*DHD);
    QRED=QRED+QADD;

    // Update STEP and HS.
    double GGSAV=GG;
    GG=ZERO;
    for(i=0;i<N;i++){
      STEP[i]=STEP[i]+ALPHA*D[i];
      HS[i]=HS[i]+ALPHA*HD[i];
      GG=GG+pow((G[i]+HS[i]),2);
    }

    // Begin another conjugate direction iteration if required.
    if(ALPHA<BSTEP){
      if(QADD<=0.01*QRED){goto n160;}
      if(GG<=1.0e-4*GGBEG){goto n160;}
      if(ITERC==ITERMAX){goto n160;}
      TEMP=GG/GGSAV;
      DD=ZERO;
      DS=ZERO;
      SS=ZERO;
      for(i=0;i<N;i++){
        D[i]=TEMP*D[i]-G[i]-HS[i];
        DD=DD+D[i]*D[i];
        DS=DS+D[i]*STEP[i];
        SS=SS+STEP[i]*STEP[i];
      }
      if(DS<=ZERO){goto n160;}
      if(SS<DELSQ){goto n40;}
    }
    CRVMIN=ZERO;
    ITERSW=ITERC;

    // Test whether an alternative iteration is required.
n90:
    if(GG<=1.0e-4*GGBEG){goto n160;}
    double SG=ZERO, SHS=ZERO, SGK;
    for(i=0;i<N;i++){
      SG=SG+STEP[i]*G[i];
      SHS=SHS+STEP[i]*HS[i];
    }
    SGK=SG+SHS;
    double ANGTEST=SGK/sqrt(GG*DELSQ);
    if(ANGTEST<=-0.99){goto n160;}

    // Begin the alternative iteration by calculating D and HD and some scalar products.
    ITERC=ITERC+1;
    TEMP=sqrt(DELSQ*GG-SGK*SGK);
    double TEMPA=DELSQ/TEMP;
    double TEMPB=SGK/TEMP;
    for(i=0;i<N;i++){D[i]=TEMPA*(G[i]+HS[i])-TEMPB*STEP[i];}
    goto n170;
n120:    
    double DG=ZERO,DHS=ZERO;
    DHD=ZERO;
    for(i=0;i<N;i++){
      DG=DG+D[i]*G[i];
      DHD=DHD+HD[i]*D[i];
      DHS=DHS+HD[i]*STEP[i];
    }

    // Seek the value of the angle that minimizes Q.
    double CF=HALF*(SHS-DHD);
    double QBEG=SG+CF;
    double QSAV=QBEG;
    double QMIN=QBEG;
    double ANGLE, CTH, STH, QNEW;
    int ISAVE=0,IU=49;
    TEMP=TWOPI/(double)(IU+1);
    for(i=1;i<=IU;i++){
      ANGLE=(double)i*TEMP;
      CTH=cos(ANGLE);
      STH=sin(ANGLE);
      QNEW=(SG+CF*CTH)*CTH+(DG+DHS*CTH)*STH;
      if(QNEW<QMIN){
        QMIN=QNEW;
        ISAVE=i;
        TEMPA=QSAV;
      }
      else if(i == ISAVE+1){
        TEMPB=QNEW;
      }
      QSAV=QNEW;
    }
    if(ISAVE==ZERO){TEMPA=QNEW;}
    if(ISAVE==IU){TEMPB=QBEG;}
    ANGLE=ZERO;
    if(TEMPA!=TEMPB){
      TEMPA=TEMPA-QMIN;
      TEMPB=TEMPB=QMIN;
      ANGLE=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB);
    }
    ANGLE=TEMP*((double)ISAVE+ANGLE);
    
    // Calculate the new STEP and HS. Then test for convergence.
    CTH=cos(ANGLE);
    STH=sin(ANGLE);
    double REDUC=QBEG-(SG+CF*CTH)*CTH-(DG+DHS*CTH)*STH;
    GG=ZERO;
    for(i=0;i<N;i++){
      STEP[i]=CTH*STEP[i]+STH*D[i];
      HS[i]=CTH*HS[i]+STH*HD[i];
      GG=GG+pow((G[i]+HS[i]),2);
    }
    QRED=QRED+REDUC;
    double RATIO=REDUC/QRED;
    if(ITERC<ITERMAX && RATIO>0.01){goto n90;}
n160:
    return 0;

    /* The following instructions act as a subroutine for setting the vector
     HD to the vector D multiplied by the second derivative matrix of Q.
     They are called from three different places, which are distinguished
     by the value of ITERC. */
n170:  
    for(i=0;i<N;i++){HD[i]=ZERO;}
    for(k=0;k<NPT;k++){
      TEMP=ZERO;
      for(j=0;j<N;j++){
        TEMP=TEMP+XPT[k][j]*D[j];
        
      }
      TEMP=TEMP*PQ[k];
      for(i=0;i<N;i++){
        HD[i]=HD[i]+TEMP*XPT[k][i];
      }
    }
    IH=0;
    for(j=0;j<=N-1;j++){
      for(i=0;i<=j;i++){
        IH=IH+1;
        if(i<j){HD[j]=HD[j]+HQ[IH-1]*D[i];}
        HD[i]=HD[i]+HQ[IH-1]*D[j];
      }
    }
    if(ITERC==0){goto n20;}
    if(ITERC<=ITERSW){goto n50;}
    goto n120;
}