#include "newuoa.h"

int newuob(double *WW, int N, int NPT, double *X, double RHOBEG, double RHOEND, int IPRINT, int MAXFUN, double *XBASE, double *XOPT, double *XNEW, double *WXPT, double *FVAL, double *GQ, double *HQ, double *PQ, double *WBMAT, double *WZMAT, int NDIM, double *D, double *VLAG, double *W){
    /*The arguments N, NPT, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical
       to the corresponding arguments in SUBROUTINE NEWUOA.
     XBASE will hold a shift of origin that should reduce the contributions
       from rounding errors to values of the model and Lagrange functions.
     XOPT will be set to the displacement from XBASE of the vector of
       variables that provides the least calculated F so far.
     XNEW will be set to the displacement from XBASE of the vector of
       variables for the current calculation of F.
     XPT will contain the interpolation point coordinates relative to XBASE.
     FVAL will hold the values of F at the interpolation points.
     GQ will hold the gradient of the quadratic model at XBASE.
     HQ will hold the explicit second derivatives of the quadratic model.
     PQ will contain the parameters of the implicit second derivatives of
       the quadratic model.
     BMAT will hold the last N columns of H.
     ZMAT will hold the factorization of the leading NPT by NPT submatrix of
       H, this factorization being ZMAT times Diag(DZ) times ZMAT^T, where
       the elements of DZ are plus or minus one, as specified by IDZ.
     NDIM is the first dimension of BMAT and has the value NPT+N.
     D is reserved for trial steps from XOPT.
     VLAG will contain the values of the Lagrange functions at a new point X.
       They are part of a product that requires VLAG to be of length NDIM.
     The array W will be used for working space. Its length must be at least
       10*NDIM = 10*(NPT+N).*/

    // Set some constants.
    double TEMPQ,SUM,RHO,DELTA;
    int IP,KNEW;
    double CRVMIN,DSQ,DNORM,RATIO;
    int ITEMP,JPT,IPT;
    double XIPT,XJPT;
    double TEMP;
    double FBEG,FOPT;
    int KOPT;
    double **XPT,**BMAT,**ZMAT,**Y;
    BMAT=matrix_reshape2d(NDIM,WW);
    XPT=matrix_reshape2d(NPT,WW);
    ZMAT=matrix_reshape2d(NPT,WW);
    Y=matrix2d(10,10);
    

    double HALF=0.5,ONE=1.0,TENTH=0.1,ZERO=0.0,F,DSTEP,DIFFC;
    int NP=N+1;
    int NH=(N*NP)/2;
    int NPTM=NPT-NP;
    int NFTEST=MAXFUN>1?MAXFUN:1;
    int i,j,k,IH;

    // Set the initial elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
    
    for(j=0;j<N;j++){
        XBASE[j]=X[j];
        for(k=0;k<NPT;k++){XPT[k][j]=ZERO;}
        for(i=0;i<NDIM;i++){BMAT[i][j]=ZERO;}
    }
    
    for(IH=0;IH<NH;IH++){HQ[IH]=ZERO;}
    for(k=0;k<NPT;k++){
        PQ[k]=ZERO;
        for(j=0;j<NPTM;j++){ZMAT[k][j]=ZERO;}
    }
    

    // Begin the initialization procedure. NF becomes one more than the number
    // of function values so far. The coordinates of the displacement of the
    // next initial interpolation point from XBASE are set in XPT(NF,.).
    double RHOSQ,RECIP,RECIQ;
    RHOSQ=RHOBEG*RHOBEG;
    RECIP=ONE/RHOSQ;
    RECIQ=sqrt(HALF)/RHOSQ;
    int NF=0;

n50:
    int NFM=NF;
    int NFMM=NF-N;
    NF=NF+1;
    if(NFM<=2*N){
        if(NFM>=1 && NFM<=N){
            XPT[NF-1][NFM-1]=RHOBEG;
        }
        else if(NFM>N){XPT[NF-1][NFMM-1]=-RHOBEG;}
    }
    else{
        ITEMP=(NFMM-1)/N;
        JPT=NFM-ITEMP*N-N;
        IPT=JPT+ITEMP;
        if(IPT>N){
            ITEMP=JPT;
            JPT=IPT-N;
            IPT=ITEMP;
        }
        XIPT=RHOBEG;
        if(FVAL[IPT+NP-1]<FVAL[IPT+1-1]){XIPT=-XIPT;}
        XJPT=RHOBEG;
        if(FVAL[JPT+NP-1]<FVAL[JPT+1-1]){XJPT=-XJPT;}
        XPT[NF-1][IPT-1]=XIPT;
        XPT[NF-1][JPT-1]=XJPT;
    }

    // Calculate the next value of F, label 70 being reached immediately
    // after this calculation. The least function value so far and its index
    // are required.
    for(j=0;j<N;j++){X[j]=XPT[NF-1][j]+XBASE[j];}
     
    goto n310;
    
n70:
    FVAL[NF-1]=F;
    if(NF==1){
        FBEG=F;
        KOPT=1;
        FOPT=F;
    }
    else if(F<FOPT){
        FOPT=F;
        KOPT=NF;
    }

    // Set the nonzero initial elements of BMAT and the quadratic model in
    // the cases when NF is at most 2*N+1.
  
    if(NFM<=2*N){
        if(NFM>=1 && NFM<=N){
            GQ[NFM-1]=(F-FBEG)/RHOBEG;
            if(NPT<NF+N){
                BMAT[0][NFM-1]=-ONE/RHOBEG;
                BMAT[NF-1][NFM-1]=ONE/RHOBEG;
                BMAT[NPT+NFM-1][NFM-1]=-HALF*RHOSQ;
            }
        }
        else if(NFM>N){
            BMAT[NF-N-1][NFMM-1]=HALF/RHOBEG;
            BMAT[NF-1][NFMM-1]=-HALF/RHOBEG;
            ZMAT[0][NFMM-1]=-RECIQ-RECIQ;
            ZMAT[NF-N-1][NFMM-1]=RECIQ;
            ZMAT[NF-1][NFMM-1]=RECIQ;
            IH=((NFMM)*(NFMM+1))/2;
            TEMP=(FBEG-F)/RHOBEG;
            HQ[IH-1]=(GQ[NFMM-1]-TEMP)/RHOBEG;
            GQ[NFMM-1]=HALF*(GQ[NFMM-1]+TEMP);
        }
    }

    // Set the off-diagonal second derivatives of the Lagrange functions and
    // the initial quadratic model.
    else{
        IH=(IPT*(IPT-1))/2+JPT;
        if(XIPT<ZERO){IPT=IPT+N;}
        if(XJPT<ZERO){JPT=JPT+N;}
        ZMAT[0][NFMM-1]=RECIP;
        ZMAT[NF-1][NFMM-1]=RECIP;
        ZMAT[IPT+1-1][NFMM-1]=-RECIP;
        ZMAT[JPT+1-1][NFMM-1]=-RECIP;
        HQ[IH-1]=(FBEG-FVAL[IPT+1-1]-FVAL[JPT+1-1]+F)/(XIPT*XJPT);
    }
    if(NF<NPT){goto n50;}

    // Begin the iterative procedure, because the initial model is complete.
    RHO=RHOBEG;
    DELTA=RHO;
    int IDZ=1;
    double DIFFA=ZERO;
    double DIFFB=ZERO;
    int ITEST=0;
    double XOPTSQ=ZERO;
    for(i=0;i<N;i++){
        XOPT[i]=XPT[KOPT-1][i];
        XOPTSQ=XOPTSQ+XOPT[i]*XOPT[i];
    }

n90:
    int NFSAV=NF;

    // Generate the next trust region step and test its length. Set KNEW
    // to -1 if the purpose of the next F will be to improve the model.
n100:
    KNEW=0;
    trsapp(N,NPT,XOPT,XPT,GQ,HQ,PQ,DELTA,D,W,slice1d(W,NP-1,10001),slice1d(W,NP+N-1,10001),slice1d(W,NP+2*N-1,10001),CRVMIN); 
    DSQ=ZERO;

    for(i=0;i<N;i++){DSQ=DSQ+D[i]*D[i];}
    DNORM=DELTA<sqrt(DSQ)?DELTA:sqrt(DSQ);
    if(DNORM<HALF*RHO){
        KNEW=-1;
        DELTA=TENTH*DELTA;
        RATIO=-1.0;
        if(DELTA<=1.5*RHO){DELTA=RHO;}
        if(NF<=NFSAV+2){goto n460;}
        TEMP=0.125*CRVMIN*RHO*RHO;
        if(TEMP<=max3d(DIFFA,DIFFB,DIFFC)){goto n460;}
       
        goto n490;
    }

    // Shift XBASE if XOPT may be too far from XBASE. First make the changes
    // to BMAT that do not depend on ZMAT.
n120:
    if(DSQ<=1.0e-3*XOPTSQ){
        TEMPQ=0.25*XOPTSQ;
        for(k=0;k<NPT;k++){
            SUM=ZERO;
            for(i=0;i<N;i++){SUM=SUM+XPT[k][i]*XOPT[i];}
            TEMP=PQ[k]*SUM;
            SUM=SUM-HALF*XOPTSQ;
            W[NPT+k]=SUM;
            for(i=0;i<=N-1;i++){
                GQ[i]=GQ[i]+TEMP*XPT[k][i];
                XPT[k][i]=XPT[k][i]-HALF*XOPT[i];
                VLAG[i]=BMAT[k][i];
                W[i]=SUM*XPT[k][i]+TEMPQ*XOPT[i];
                IP=NPT+i;
                for(j=0;j<=i;j++){
                    BMAT[IP][j]=BMAT[IP][j]+VLAG[i]*W[j]+W[i]*VLAG[j];
                }
            }
        }

        // Then the revisions of BMAT that depend on ZMAT are calculated.
        double SUMZ;
        for(k=0;k<NPTM;k++){
            SUMZ=ZERO;
            for(i=0;i<NPT;i++){
                SUMZ=SUMZ+ZMAT[i][k];
                W[i]=W[NPT+i]*ZMAT[i][k];
            }
            for(j=0;j<N;j++){
                SUM=TEMPQ*SUMZ*XOPT[j];
                for(i=0;i<NPT;i++){SUM=SUM+W[i]*XPT[i][j];}
                VLAG[j]=SUM;
                if(k<IDZ-1){SUM=-SUM;}
                for(i=0;i<NPT;i++){BMAT[i][j]=BMAT[i][j]+SUM*ZMAT[i][k];}
            }
            for(i=0;i<=N-1;i++){
                IP=i+NPT;
                TEMP=VLAG[i];
                if(k<IDZ-1){TEMP=-TEMP;}
                for(j=0;j<=i;j++){BMAT[IP][j]=BMAT[IP][j]+TEMP*VLAG[j];}
            }
        }

        // The following instructions complete the shift of XBASE, including
        // the changes to the parameters of the quadratic model.
        IH=0;
        for(j=0;j<=N-1;j++){
            W[j]=ZERO;
            for(k=0;k<NPT;k++){
                W[j]=W[j]+PQ[k]*XPT[k][j];
                XPT[k][j]=XPT[k][j]-HALF*XOPT[j];
            }
            for(i=0;i<=j;i++){
                IH=IH+1;
                if(i<j){GQ[j]=GQ[j]+HQ[IH-1]*XOPT[i];}
                GQ[i]=GQ[i]+HQ[IH-1]*XOPT[j];
                HQ[IH-1]=HQ[IH-1]+W[i]*XOPT[j]+XOPT[i]*W[j];
                BMAT[NPT+i][j]=BMAT[NPT+j][i];
            }
        }
        for(j=0;j<N;j++){
            XBASE[j]=XBASE[j]+XOPT[j];
            XOPT[j]=ZERO;
        }
        XOPTSQ=ZERO;
    }

    // Pick the model step if KNEW is positive. A different choice of D
    // may be made later, if the choice of D by BIGLAG causes substantial
    // cancellation in DENOM.
    double ALPHA;
    
    if(KNEW>0){
        ALPHA=biglag(N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ-1,NDIM,KNEW-1,DSTEP,D,ALPHA,VLAG,slice1d(VLAG,NPT+1-1,10001),W,slice1d(W,NP-1,10001),slice1d(W,NP+N-1,10001));
    }
    

    // Calculate VLAG and BETA for the current choice of D. The first NPT
    // components of W_check will be held in W.
    double SUMA,SUMB;
    for(k=0;k<NPT;k++){
        SUMA=ZERO;
        SUMB=ZERO;
        SUM=ZERO;
        for(j=0;j<N;j++){
            SUMA=SUMA+XPT[k][j]*D[j];
            SUMB=SUMB+XPT[k][j]*XOPT[j];
            SUM=SUM+BMAT[k][j]*D[j];
        }
        W[k]=SUMA*(HALF*SUMA+SUMB);
        VLAG[k]=SUM;
    }
    double BETA=ZERO;
    for(k=0;k<NPTM;k++){
        SUM=ZERO;
        for(i=0;i<NPT;i++){SUM=SUM+ZMAT[i][k]*W[i];}
        if(k<IDZ-1){
            BETA=BETA+SUM*SUM;
            SUM=-SUM;
        }
        else{BETA=BETA-SUM*SUM;}
        
        for(i=0;i<NPT;i++){VLAG[i]=VLAG[i]+SUM*ZMAT[i][k];}
    }


    double BSUM=ZERO;
    double DX=ZERO;
    int JP;
    for(j=0;j<N;j++){
       
        SUM=ZERO;
        for(i=0;i<NPT;i++){SUM=SUM+W[i]*BMAT[i][j];}
        BSUM=BSUM+SUM*D[j];
        JP=NPT+j;
        for(k=0;k<N;k++){SUM=SUM+BMAT[JP][k]*D[k];}
        VLAG[JP]=SUM;
        BSUM=BSUM+SUM*D[j];
        DX=DX+D[j]*XOPT[j];
    }
    BETA=DX*DX+DSQ*(XOPTSQ+DX+DX+HALF*DSQ)+BETA-BSUM;
    VLAG[KOPT-1]=VLAG[KOPT-1]+ONE;

    // If KNEW is positive and if the cancellation in DENOM is unacceptable,
    // then BIGDEN calculates an alternative model step, XNEW being used for
    // working space.
    if(KNEW>0){
        TEMP=ONE+ALPHA*BETA/(VLAG[KNEW-1]*VLAG[KNEW-1]);

        if(fabs(TEMP)<=0.8){
            bigden(N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ-1,NDIM,KOPT-1,KNEW-1,D,W,VLAG,BETA,XNEW,NDIM+1,6*NDIM+1);
        }
        
    }

    // Calculate the next value of the objective function.
n290:
    for(i=0;i<N;i++){
        XNEW[i]=XOPT[i]+D[i];
        X[i]=XBASE[i]+XNEW[i];
    }
    
    
    NF=NF+1;
n310:

    if(NF>NFTEST){
        NF=NF-1;
        if(IPRINT>0){
            printf("Return from NEWUOA because CALFUN has been called MAXFUN times.\n");
        }
        goto n530;
    }

    F=calfun(N,X,F,Y);

    if(IPRINT==3){
        printf("Function number: %d, F=%lf, The corresponding X is:\n",NF,F);
        for(i=0;i<N;i++){
            printf("%lf ",X[i]);
        }
        printf("\n");
    }
    
    if(NF<=NPT){goto n70;}
    if(KNEW==-1){goto n530;}

    // Use the quadratic model to predict the change in F due to the step D,
    // and set DIFF to the error of this prediction.
    double VQUAD=ZERO;
    IH=0;
    for(j=0;j<=N-1;j++){
        VQUAD=VQUAD+D[j]*GQ[j];
        for(i=0;i<=j;i++){
            IH=IH+1;
            TEMP=D[i]*XNEW[j]+D[j]*XOPT[i];
            if(i==j){TEMP=HALF*TEMP;}
            VQUAD=VQUAD+TEMP*HQ[IH-1];
        }
    }
    for(k=0;k<NPT;k++){VQUAD=VQUAD+PQ[k]*W[k];}
    double DIFF;
    DIFF=F-FOPT-VQUAD;
    DIFFC=DIFFB;
    DIFFB=DIFFA;
    DIFFA=fabs(DIFF);
    if(DNORM>RHO){NFSAV=NF;}

    // Update FOPT and XOPT if the new F is the least value of the objective
    // function so far. The branch when KNEW is positive occurs if D is not
    // a trust region step.
    double FSAVE=FOPT;
    if(F<FOPT){
        FOPT=F;
        XOPTSQ=ZERO;
        for(i=0;i<N;i++){
            XOPT[i]=XNEW[i];
            XOPTSQ=XOPTSQ+XOPT[i]*XOPT[i];
        }
    }
    int KSAVE=KNEW;
    if(KNEW>0){goto n410;}

    // Pick the next value of DELTA after a trust region step.
    if(VQUAD>=ZERO){
        if(IPRINT>0){
            printf("Return from NEWUOA because a trust region step has failed to reduce Q.\n");
        }
        goto n530;
    }
    RATIO=(F-FSAVE)/VQUAD;
    if(RATIO<=TENTH){DELTA=HALF*DNORM;}
    else if(RATIO<=0.7){DELTA=HALF*DELTA>DNORM?HALF*DELTA:DNORM;}
    else{DELTA=HALF*DELTA>DNORM*2?HALF*DELTA:DNORM*2;}
    if(DELTA<=1.5*RHO){DELTA=RHO;}


    //Set KNEW to the index of the next interpolation point to be deleted.
    RHOSQ=pow(max2d(TENTH*DELTA,RHO),2);
    int KTEMP=0;
    double DETRAT=ZERO,HDIAG,DISTSQ;
    if(F>=FSAVE){
        KTEMP=KOPT;
        DETRAT=ONE;
    }
    for(k=0;k<NPT;k++){
        HDIAG=ZERO;
        for(j=0;j<NPTM;j++){
            TEMP=ONE;
            if(j<IDZ-1){TEMP=-ONE;}
            HDIAG=HDIAG+TEMP*ZMAT[k][j]*ZMAT[k][j];
        }
        TEMP=fabs(BETA*HDIAG+VLAG[k]*VLAG[k]);
        DISTSQ=ZERO;
        for(j=0;j<N;j++){DISTSQ=DISTSQ+pow((XPT[k][j]-XOPT[j]),2);}
        if(DISTSQ>RHOSQ){TEMP=TEMP*pow((DISTSQ/RHOSQ),3);}
        if(TEMP>DETRAT && k+1!=KTEMP){
            DETRAT=TEMP;
            KNEW=k+1;
        }
        continue;
    }
    if(KNEW==0){goto n460;}

    // Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
    // can be moved. Begin the updating of the quadratic model, starting
    // with the explicit second derivative term.
n410:
    update(N,NPT,BMAT,ZMAT,IDZ-1,NDIM,VLAG,BETA,KNEW-1,W);

    FVAL[KNEW-1]=F;
    IH=0;
    for(i=0;i<=N-1;i++){
        TEMP=PQ[KNEW-1]*XPT[KNEW-1][i];
        for(j=0;j<=i;j++){
            IH=IH+1;
            HQ[IH-1]=HQ[IH-1]+TEMP*XPT[KNEW-1][j];
        }
    }
    PQ[KNEW-1]=ZERO;

    // Update the other second derivative parameters, and then the gradient
    // vector of the model. Also include the new interpolation point.
    for(j=0;j<NPTM;j++){
        TEMP=DIFF*ZMAT[KNEW-1][j];
        if(j<IDZ-1){TEMP=-TEMP;}
        for(k=0;k<NPT;k++){PQ[k]=PQ[k]+TEMP*ZMAT[k][j];}
    }
    double GQSQ=ZERO;
    for(i=0;i<N;i++){
        GQ[i]=GQ[i]+DIFF*BMAT[KNEW-1][i];
        GQSQ=GQSQ+GQ[i]*GQ[i];
        XPT[KNEW-1][i]=XNEW[i];
    }

    // If a trust region step makes a small change to the objective function,
    // then calculate the gradient of the least Frobenius norm interpolant at
    // XBASE, and store it in W, using VLAG for a vector of right hand sides.
    double GISQ;
    if(KSAVE==0 && DELTA==RHO){
        if(fabs(RATIO)>1.0e-2){
            ITEST=0;
        }
        else{
            for(i=0;i<NPT;i++){VLAG[k]=FVAL[k]-FVAL[KOPT-1];}
            GISQ=ZERO;
            for(i=0;i<N;i++){
                SUM=ZERO;
                for(k=0;k<NPT;k++){SUM=SUM+BMAT[k][i]*VLAG[k];}
                GISQ=GISQ+SUM*SUM;
                W[i]=SUM;
            }
            // Test whether to replace the new quadratic model by the least Frobenius
            // norm interpolant, making the replacement if the test is satisfied.
            ITEST=ITEST+1;
            if(GQSQ<1.0e2*GISQ){ITEST=0;}
            if(ITEST>=3){
                for(i=0;i<N;i++){GQ[i]=W[i];}
                for(IH=0;IH<NH;IH++){HQ[IH]=ZERO;}
                for(j=0;j<NPTM;j++){
                    W[j]=ZERO;
                    for(k=0;k<NPT;k++){
                        W[j]=W[j]+VLAG[k]*ZMAT[k][j];
                    }
                    if(j<IDZ-1){W[j]=-W[j];}
                }
                for(k=0;k<NPT;k++){
                    PQ[k]=ZERO;
                    for(j=0;j<NPTM;j++){
                        PQ[k]=PQ[k]+ZMAT[k][j]*W[j];
                        
                    }
                }
                ITEST=0;
            }
        }

    }
    if(F<FSAVE){KOPT=KNEW;}

    // If a trust region step has provided a sufficient decrease in F, then
    // branch for another trust region calculation. The case KSAVE>0 occurs
    // when the new function value was calculated by a model step.
    if(F<=FSAVE+TENTH*VQUAD){goto n100;}
    if(KSAVE>0){goto n100;}

    // Alternatively, find out if the interpolation points are close enough
    // to the best point so far.
    KNEW=0;
n460:
    DISTSQ=4.0*DELTA*DELTA;
    for(k=0;k<NPT;k++){
        SUM=ZERO;
        for(j=0;j<N;j++){SUM=SUM+pow((XPT[k][j]-XOPT[j]),2);}
        if(SUM>DISTSQ){
            KNEW=k+1;
            DISTSQ=SUM;
        }
        continue;
    }

    // If KNEW is positive, then set DSTEP, and branch back for the next
    // iteration, which will generate a "model step".
    if(KNEW>0){
        DSTEP=max2d(min2d(TENTH*sqrt(DISTSQ),HALF*DELTA),RHO);
        DSQ=DSTEP*DSTEP;
        goto n120;
    }
    if(RATIO>ZERO){goto n100;}
    if(max2d(DELTA,DNORM)>RHO){goto n100;}

    // The calculations with the current value of RHO are complete. Pick the
    // next values of RHO and DELTA.
n490:
    if(RHO>RHOEND){
        DELTA=HALF*RHO;
        RATIO=RHO/RHOEND;
        if(RATIO<=16.0){RHO=RHOEND;}
        else if(RATIO<=250.0){RHO=sqrt(RATIO)*RHOEND;}
        else{RHO=TENTH*RHO;}
        DELTA=max2d(DELTA,RHO);
        if(IPRINT>=2){
            if(IPRINT>=3){printf("     ");}
            printf("New RHO=%lf,Number of function values=%d\n",RHO,NF);
            printf("Least value of F=%lf, The corresponding X is:\n",FOPT);
            for(i=0;i<N;i++){
            printf("%lf ",XBASE[i]+XOPT[i]);
            }
            printf("\n");
        }
        goto n90;
    }


// Return from the calculation, after another Newton-Raphson step, if
    // it is too short to have been tried before.
    if(KNEW==-1){goto n290;}
n530:
    if(FOPT<=F){
        for(i=0;i<N;i++){X[i]=XBASE[i]+XOPT[i];}
        F=FOPT;
    }
    if(IPRINT>=1){
        printf("  \n");
        printf("At the return from NEWUOA  Number of function value = %d\n",NF);
        printf("Least value of F=%lf, The corresponding X is:\n",F);
        for(i=0;i<N;i++){
            printf("%lf ",X[i]);
            }
            printf("\n");
        printf("       \n\n");
    }

    return 0;
}