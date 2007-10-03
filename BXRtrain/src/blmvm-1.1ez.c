/* 
COPYRIGHT NOTIFICATION

(C) COPYRIGHT 2004 UNIVERSITY OF CHICAGO

This program discloses material protectable under copyright laws of the United States.
Permission to copy and modify this software and its documentation is
hereby granted, provided that this notice is retained thereon and on all copies or
modifications. The University of Chicago makes no representations as to the suitability
and operability of this software for any purpose. It is provided "as is"; without
express or implied warranty. Permission is hereby granted to use, reproduce, prepare
derivative works, and to redistribute to others, so long as this original copyright notice
is retained.  

Any publication resulting from research that made use of this software 
should cite the document:  Steven J. Benson and Jorge Mor\'{e}, 
"A Limited-Memory Variable-Metric Algorithm for Bound-Constrained Minimization",
Mathematics and Computer Science Division, Argonne National Laboratory,
ANL/MCS-P909-0901, 2001.




    Argonne National Laboratory with facilities in the states of Illinois and Idaho, is
    owned by The United States Government, and operated by the University of Chicago under
    provision of a contract with the Department of Energy.

                                    DISCLAIMER
    THIS PROGRAM WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY AN AGENCY OF THE UNITED
    STATES GOVERNMENT. NEITHER THE UNITED STATES GOVERNMENT NOR ANY AGENCY THEREOF, NOR THE
    UNIVERSITY OF CHICAGO, NOR ANY OF THEIR EMPLOYEES OR OFFICERS, MAKES ANY WARRANTY, EXPRESS
    OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
    COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED,
    OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. REFERENCE HEREIN TO
    ANY SPECIFIC COMMERCIAL PRODUCT, PROCESS, OR SERVICE BY TRADE NAME, TRADEMARK,
    MANUFACTURER, OR OTHERWISE, DOES NOT NECESSARILY CONSTITUTE OR IMPLY ITS ENDORSEMENT,
    RECOMMENDATION, OR FAVORING BY THE UNITED STATES GOVERNMENT OR ANY AGENCY THEREOF. THE
    VIEW AND OPINIONS OF AUTHORS EXPRESSED HEREIN DO NOT NECESSARILY STATE OR REFLECT THOSE OF
    THE UNITED STATES GOVERNMENT OR ANY AGENCY THEREOF.
*/ 
/*
   Applications of the BLMVM solver for bound constrained minimization must 
   implement 2 routines: BLMVMFunctionGradient(), BLMVMConverge(). In addition, they must
   call the routine BLMVMSolveIt() with the number of variables, and initial solution,
   lower and upper bounds, and a parameter.  To solve applications other than the following example,
   replace the two routines with other routines and call BLMVMSolveIt() with the appropriate arguments.
*/
#include <math.h> 
#include <stdlib.h> 
#include <stdio.h> 
#include <string.h> 
#ifndef BLMVMABC_H
#define BLMVMABC_H
extern int BLMVMFunctionGradient(void*,double[],double[],int,double*,double*);
extern int BLMVMConverge(void*,int,double,double[],int);
int BLMVMSolveIt(void*,double[],double[],double[],int,int);
static int BLMVMSolve(int(*)(void*,double[],double[],int,double*,double*),
		      int (*)(void*,int,double,double[],int),
		      void*,double[],double[],double[],int,int);
#endif
/* ------------ Begin blmvmsolveit.c ----------------- */

#undef __FUNCT__  
#define __FUNCT__ "BLMVMSolveIt"
/* --------------------------------------------------------------------------
   BLMVMSolveIt -- User must implement the routine that evaluates 
   the objective function and its gradient at current point x 

   Input:   ctx - pointer to user-defined context.
            xl - lower bounds on variables (array of length n).
            x - current variables (array of length n).
            xu - upper bounds variables (array of length n).
            n - number of variables.
            lm - number of correction pairs to be used in the solver (>0)(Suggested value: 5).
            maxstep - maximum change in solution per iteration (>0)(Suggested value: 1e10)
            f - address of current objective function value;

   Output:  ctx - unchanged
            xl - sensitivity of objective to the lower bounds (array of length n).
            x - solution (array of length n).
            xu - sensitivity of objective to the upper bounds (array of length n).

   return:  0 for Normal return and nonzero to abort mission.
*/
int BLMVMSolveIt(void* ctx, double xl[], double x[], double xu[], int n, int lm){
  int info; 
  /* Set solver parameters and solve. */
  /* Creation of app struct */
  info=BLMVMSolve(BLMVMFunctionGradient,BLMVMConverge,ctx,xl,x,xu,n,lm);
  return(info);
}

/* ------------ End blmvmsolveit.c ----------------- */
/* --------- Code below this line does not have to be modified -----------------*/


#ifndef BLMVMERR_H
#define BLMVMERR_H

//#define __FUNCT__ "Function Unknown"

#ifndef CHKERRQ
#define CHKERRQ(a)  { if (a){ printf("Error detected in "); printf(__FUNCT__); return a; } }
#endif

#define BLMVMCALLOC(VAR,SIZE,TYPE)  {VAR = (TYPE*)calloc(SIZE, sizeof(TYPE));}
#define BLMVMFREE(VAR)              {if (VAR){free(VAR);}}

#endif
/* Begin blmvm.h -- define Vector interface used by BLMVM solver, matrix, and line search */
#ifndef BLMVMVECTOR_H
#define BLMVMVECTOR_H

typedef struct _P_BLMVMVec *BLMVMVec;

static int BLMVM_StepBound(BLMVMVec, BLMVMVec, BLMVMVec,BLMVMVec, double*,double*,double*);
static int BLMVMVecDuplicate(BLMVMVec,BLMVMVec*);
static int BLMVMVecCopy(BLMVMVec,BLMVMVec);
static int BLMVMVecAXPY(double,BLMVMVec,BLMVMVec);
static int BLMVMVecMultiply(double,BLMVMVec,BLMVMVec, double, double*, double*);
static int BLMVMVecAYPX(double,BLMVMVec,BLMVMVec);
static int BLMVMVecWAXPY(double,BLMVMVec,BLMVMVec,BLMVMVec);
static int BLMVMVecApplyBounds(BLMVMVec,BLMVMVec,BLMVMVec);
static int BLMVMVecDot(BLMVMVec,BLMVMVec,double*);
static int BLMVMVec2Norm(BLMVMVec,double*);
static int BLMVMVecDestroy(BLMVMVec*);
static int BLMVMVecScale(double,BLMVMVec);
static int BLMVMVecProjectGradient(BLMVMVec,BLMVMVec,BLMVMVec,BLMVMVec,BLMVMVec,BLMVMVec,BLMVMVec);
static int BLMVMVecCompatible(BLMVMVec,BLMVMVec, int*);
static int BLMVMVecView(BLMVMVec V);
static int BLMVMVecDim(BLMVMVec V);
#endif
/* End blmvm.h */
/* --------- Begin blmvmapplication.h ---------- */

#ifndef BLMVMAPPLICATION_H
#define BLMVMAPPLICATION_H

typedef struct _P_BLMVM_APPLICATION* BLMVM_APP;
static int BLMVMComputeFunctionGradient(BLMVM_APP, BLMVMVec,double*,BLMVMVec,BLMVMVec);
static int BLMVMComputeBounds(BLMVM_APP, BLMVMVec, BLMVMVec);
static int BLMVMMonitor(BLMVM_APP,int,double,double,BLMVMVec);

#endif
/* --------- End blmvmapplication.h ---------- */
/* -----    Begin blmvm.h ------------- */
#ifndef BLMVMSOLVER_H
#define BLMVMSOLVER_H


typedef struct  P_BLMVM *BLMVM;

static int BLMVMCreate(BLMVM*);
static int BLMVMSetUp(BLMVM, BLMVMVec);
static int BLMVMSetHistory(BLMVM,int);
static int BLMVMGetHistory(BLMVM,int*);
static int BLMVMMinimize(BLMVM, BLMVM_APP);
static int BLMVMGetProjectedGradientVec(BLMVM, BLMVMVec*);
static int BLMVMGetGradientVec(BLMVM, BLMVMVec*);
static int BLMVMGetBoundDualVec(BLMVM, BLMVMVec*, BLMVMVec*);
static int BLMVMGetStepVec(BLMVM, BLMVMVec*);
static int BLMVMGetSolutionVec(BLMVM,BLMVMVec*);
static int BLMVMGetStopReason(BLMVM,int*);
static int BLMVMTakeDown(BLMVM);
static int BLMVMDestroy(BLMVM*);
static int BLMVMSetMaxStepNorm(BLMVM,double);
static int BLMVMGetStepNorm(BLMVM,double*);
static int BLMVMGetStepSize(BLMVM,double*);
#endif
/* -----    End blmvm.h ------------- */
/* --------- Begin blmvmapplication2.h ---------- */
#ifndef BLMVMAPPLICATION2_H
#define BLMVMAPPLICATION2_H
static int BLMVMSetFunctionGradient(BLMVM_APP, int(*)(void*, double*,double*,int,double*,double*),void*);
static int BLMVMSetBounds(BLMVM_APP, int(*)(void*, double*, double*, int), void*);
static int BLMVMSetConverge(BLMVM_APP, int (*)(void*,int,double,double*,int), void*);
static int BLMVMAppCreate(BLMVM_APP *);
static int BLMVMAppDestroy(BLMVM_APP *);

static int BLMVMVecCreateWArray(BLMVMVec *VV, double*, int);
static int BLMVMVecCreate(int,BLMVMVec*);
static int BLMVMVecGetDoubles(BLMVMVec,double**,int*);

#endif
/* --------- End blmvmapplication2.h ---------- */
/* ------------ Begin blmvmabc.c ----------------- */

static int BLMVMLowerAndUpperBounds(void* ctx, double xl[], double xu[], int n);
typedef struct{
  double *xl,*xu;
  int n;
} ABCCtx;

#undef __FUNCT__  
#define __FUNCT__ "BLMVMSolve"
/* --------------------------------------------------------------------------
   BLMVMSolve -- User must implement this routine that evaluates 
   the objective function and its gradient at current point x 

   Input:   ctx - pointer to user-defined context.
            xl - lower bounds on variables (array of length n).
            x - current variables (array of length n).
            xu - upper bounds variables (array of length n).
            n - number of variables.
            lm - number of correction pairs to be used in the solver (>0)(Suggested value: 5).
            maxstep - maximum change in solution per iteration (>0)(Suggested value: 1e10)
            f - address of current objective function value;

   Output:  ctx - unchanged
            xl - sensitivity of objective to the lower bounds (array of length n).
            x - solution (array of length n).
            xu - sensitivity of objective to the upper bounds (array of length n).

   return:  0 for Normal return and nonzero to abort mission.
*/
static int BLMVMSolve(int(*fg)(void*,double[],double[],int,double*,double*),int (*conv)(void*,int,double,double[],int),void* ctx, double xl[], double x[], double xu[], int n, int lm){
  int info,i,nn; 
  double dd,maxstep=2.0;  // !!
  double *dxl,*xx,*dxu;
  BLMVM blmvm;
  BLMVM_APP blmvmapp;
  BLMVMVec X,T,DXU,DXL;
  ABCCtx bctx;

  for (i=0;i<n;i++){dd=xl[i];} for (i=0;i<n;i++){dd=x[i];}  for (i=0;i<n;i++){dd=xu[i];}  /* Check for valid data */
  for (i=0;i<n;i++){if (xl[i]>xu[i]){ printf("ERROR: Lower bound %d (%4.4e) is greater that upper bound (%4.4e)",i,xl[i],xu[i]);}}

  /* Set solver parameters and solve. */
  /* Creation of app struct */
  info= BLMVMAppCreate(&blmvmapp);CHKERRQ(info);
  bctx.xl=xl;bctx.xu=xu;bctx.n=n;
  info= BLMVMSetBounds(blmvmapp, BLMVMLowerAndUpperBounds, (void*)&bctx); CHKERRQ(info);
  info= BLMVMSetFunctionGradient(blmvmapp,fg,ctx); CHKERRQ(info);
  info= BLMVMSetConverge(blmvmapp,conv, ctx); CHKERRQ(info);
    /* Creation of blmvm struct */ 
  info=BLMVMCreate(&blmvm); CHKERRQ(info);
  info=BLMVMVecCreate(n,&X); CHKERRQ(info);
  info=BLMVMVecGetDoubles(X,&xx,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){xx[i]=x[i];}
  info=BLMVMSetHistory(blmvm,lm);CHKERRQ(info);
  info=BLMVMSetMaxStepNorm(blmvm,maxstep);CHKERRQ(info);
  info=BLMVMSetUp(blmvm,X); CHKERRQ(info);
  info=BLMVMGetGradientVec(blmvm,&T);CHKERRQ(info);
  info=BLMVMGetProjectedGradientVec(blmvm,&T);CHKERRQ(info);
  info=BLMVMGetStepVec(blmvm,&T);CHKERRQ(info);
  info=BLMVMGetSolutionVec(blmvm,&T);CHKERRQ(info);
  if (info==-1){info=BLMVMVecView(T);CHKERRQ(info);}
  info=BLMVMMinimize(blmvm, blmvmapp); CHKERRQ(info);  
  info=BLMVMGetStepNorm(blmvm,&dd); CHKERRQ(info);  
  info=BLMVMGetStopReason(blmvm,&i); CHKERRQ(info);  
  info=BLMVMGetHistory(blmvm,&i); CHKERRQ(info);  
  info=BLMVMGetStepSize(blmvm,&dd); CHKERRQ(info);
  info=BLMVMGetBoundDualVec(blmvm,&DXL, &DXU);CHKERRQ(info);
  info=BLMVMVecGetDoubles(DXL,&dxl,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){xl[i]=dxl[i];}  
  info=BLMVMVecGetDoubles(DXU,&dxu,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){xu[i]=dxu[i];}
  info=BLMVMVecGetDoubles(X,&xx,&nn);CHKERRQ(info);
  for (i=0;i<nn;i++){x[i]=xx[i];}
  info=BLMVMTakeDown(blmvm); CHKERRQ(info);
  info=BLMVMDestroy(&blmvm); CHKERRQ(info);  
  info=BLMVMAppDestroy(&blmvmapp);CHKERRQ(info);
  info=BLMVMVecDestroy(&X);CHKERRQ(info);
  return(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BLMVMLowerAndUpperBounds"
static int BLMVMLowerAndUpperBounds(void* ctx, double xl[], double xu[], int n){
  ABCCtx *bctx=(ABCCtx*)ctx;
  int i;
  for (i=0;i<n;i++){xl[i]=bctx->xl[i];}
  for (i=0;i<n;i++){xu[i]=bctx->xu[i];}
  return(0);
}
/* ------------ End blmvmabc.c ----------------- */
/* ---------- Begin blmvmlinesearch.h --------------------------- */
/* Define the interface to line searches used by the BLMVM solver */
#ifndef BLMVMLINESEARCH_H
#define BLMVMLINESEARCH_H


typedef struct _P_BLMVMLineSearch *BLMVM_LINESEARCH;

static int BLMVMApplyLineSearch(BLMVM_LINESEARCH,BLMVM_APP,BLMVMVec[9],double*,double*);
static int BLMVMLineSearchCreate(BLMVM_LINESEARCH*);
static int BLMVMLineSearchDestroy(BLMVM_LINESEARCH*);
static int BLMVMLineSearchSetUp(BLMVM_LINESEARCH,BLMVMVec);
static int BLMVMLineSearchTakeDown(BLMVM_LINESEARCH);

#endif
/* End blmvmlinesearch.h */
/* --------Begin blmvmmatrix.h ---------------------- */
/* Define Matrix interface used by BLMVM solver ------*/
#ifndef BLMVMMATRIX_H
#define BLMVMMATRIX_H


typedef struct  _P_LMVMMat *LMVMMat;

int LMVMMatCreate(int,BLMVMVec,LMVMMat*);
int LMVMMatUpdate(LMVMMat,BLMVMVec,BLMVMVec);
int LMVMMatSolve(LMVMMat,BLMVMVec,BLMVMVec,BLMVMVec);
int LMVMMatDestroy(LMVMMat);
int LMVMMatRefresh(LMVMMat);
int LMVMMatX0(LMVMMat,BLMVMVec);

#endif
/* ---------- End blmvmmatrix.h ---------------------*/
/* ----------------- Begin blmvm.c --------------- */

#ifdef BKEY
#undef BKEY
#endif
#define BKEY 2323

static int BLMVMInitialize(BLMVM);
/* THIS FILE CONTAINS A C IMPLEMENTATION OF THE BLMVM SOLVER */

/* CODE BELOW THIS POINT DEFINES THE SOLVER */

struct  P_BLMVM{

  int lm;
  int setup;
  LMVMMat  M;

  BLMVMVec W;
  BLMVMVec X;
  BLMVMVec DX;
  BLMVMVec GP;
  BLMVMVec G;
  BLMVMVec XL;
  BLMVMVec XU;
  BLMVMVec DXL;
  BLMVMVec DXU;
  BLMVMVec IHI;

  BLMVM_APP Application;
  BLMVM_LINESEARCH Line;
  int pgits;
  int stopreason;
  int key;
  double maxdx;
  double alpha;
  double dxnorm;
};

#ifdef CHKSOLVER
#undef CHKSOLVER
#endif

#define CHKSOLVER(a) {if (!(a)||((a)->key!=BKEY)){ printf("BLMVMERROR: Invalid BLMVM object\n  %s\n",__FUNCT__); return (1);}}

#undef __FUNCT__  
#define __FUNCT__ "BLMVMMinimize"
static int BLMVMMinimize(BLMVM blmvm, BLMVM_APP blmvmapp){

  BLMVMVec  W,X,G,GP;
  BLMVMVec  DX,XL,XU,DXL,DXU,IHI;
  BLMVMVec  VArray[9];
  LMVMMat   M;
  int       info,flag,iter=0;
  double    f,alpha=0,gdx,dxnorm;

  CHKSOLVER(blmvm);
  blmvm->Application=blmvmapp;
  info=BLMVMGetSolutionVec(blmvm,&X); CHKERRQ(info);
  info=BLMVMSetUp(blmvm,X); CHKERRQ(info);

  W=blmvm->W;G=blmvm->G;GP=blmvm->GP;
  DX=blmvm->DX; XL=blmvm->XL,XU=blmvm->XU;DXL=blmvm->DXL,DXU=blmvm->DXU;IHI=blmvm->IHI;
  M=blmvm->M;
  VArray[0]=XL;VArray[1]=XU;VArray[2]=DX;VArray[3]=X;
  VArray[4]=G; VArray[5]=W; VArray[6]=DXL;VArray[7]=DXU;VArray[8]=IHI;
  info=BLMVMComputeBounds(blmvm->Application,XL,XU); CHKERRQ(info);
  info=BLMVMVecApplyBounds(X,XL,XU); CHKERRQ(info);
  info=BLMVMComputeFunctionGradient(blmvm->Application,X,&f,G,IHI); CHKERRQ(info);
  blmvm->pgits=0; blmvm->stopreason=0; blmvm->dxnorm=0; blmvm->alpha=0;
 
  while (1){
    info=BLMVMVecProjectGradient(G,XL,X,XU,DXL,GP,DXU); CHKERRQ(info);

    flag=0; /* STOPPING CRITERIA */
    flag=BLMVMMonitor(blmvm->Application,iter,blmvm->alpha,blmvm->dxnorm,GP);
    if (flag) break;
    iter++;

    info=LMVMMatUpdate(M,X,GP); CHKERRQ(info);
    info=LMVMMatSolve(M,G,DX,IHI); CHKERRQ(info);
    info=BLMVMVecProjectGradient(DX,XL,X,XU,DXL,W,DXU); CHKERRQ(info);
    info=BLMVMVec2Norm(W,&dxnorm);CHKERRQ(info);
	//printf("dxnorm is %f", dxnorm);

    info=BLMVMVecDot(G,W,&gdx);CHKERRQ(info);
    if (gdx<=0){
      info=BLMVMVecCopy(G,DX);CHKERRQ(info); blmvm->pgits++;
      info=BLMVMVec2Norm(GP,&dxnorm);CHKERRQ(info);
    //	printf("Xdxnorm is %f", dxnorm);
    }

    alpha=-1.0;
    if (1.0*dxnorm>blmvm->maxdx) alpha=-(blmvm->maxdx) / dxnorm;
	printf("\nalpha is %f, limit on alpha is %f\n", alpha, -(blmvm->maxdx) / dxnorm);
    info=BLMVMApplyLineSearch(blmvm->Line,blmvm->Application,VArray,&f,&alpha);CHKERRQ(info);
    blmvm->dxnorm=fabs(dxnorm*alpha);
    blmvm->alpha=fabs(alpha);
  }
  blmvm->stopreason=flag;
  return(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BLMVMCreate"
static int BLMVMCreate(BLMVM* pblmvm){
  int info;
  BLMVM blmvm;
  BLMVMCALLOC(blmvm,1,struct P_BLMVM); 
  blmvm->key=BKEY;
  info=BLMVMInitialize(blmvm); CHKERRQ(info);
  blmvm->lm=8;
  info=BLMVMLineSearchCreate(&blmvm->Line);CHKERRQ(info);
  *pblmvm=blmvm;
  return(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BLMVMDestroy"
static int BLMVMDestroy(BLMVM* pblmvm){
  int info;
  BLMVM blmvm=0;
  if (pblmvm){blmvm=*pblmvm;}
  if (blmvm){
    CHKSOLVER(blmvm);
    info=BLMVMTakeDown(blmvm);CHKERRQ(info);
    info=BLMVMLineSearchDestroy(&blmvm->Line);CHKERRQ(info);
    blmvm->key=0;
    BLMVMFREE(blmvm);
    *pblmvm=0;
  }
  return(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BLMVMSetUp"
static int BLMVMSetUp(BLMVM blmvm, BLMVMVec X){
  int info;
  int flag;

  CHKSOLVER(blmvm);
  blmvm->X=X;
  info=BLMVMLineSearchSetUp(blmvm->Line,X);CHKERRQ(info);
  if (blmvm->setup){
    info=BLMVMVecCompatible(X,blmvm->DX,&flag);CHKERRQ(info);
    if (flag){return(0);}
  }
  info=BLMVMVecDuplicate(X,&blmvm->W); CHKERRQ(info);
  info=BLMVMVecDuplicate(X,&blmvm->DX); CHKERRQ(info);
  info=BLMVMVecDuplicate(X,&blmvm->GP); CHKERRQ(info);
  info=BLMVMVecDuplicate(X,&blmvm->G); CHKERRQ(info);
  info=BLMVMVecDuplicate(X,&blmvm->XL); CHKERRQ(info);
  info=BLMVMVecDuplicate(X,&blmvm->XU); CHKERRQ(info);  
  info=BLMVMVecDuplicate(X,&blmvm->DXL); CHKERRQ(info);
  info=BLMVMVecDuplicate(X,&blmvm->DXU); CHKERRQ(info);  
  info=BLMVMVecDuplicate(X,&blmvm->IHI); CHKERRQ(info);  
  info=LMVMMatCreate(blmvm->lm,X,&blmvm->M);CHKERRQ(info);
  blmvm->setup=1;
  return(0);
}

/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "BLMVMInitialize"
static int BLMVMInitialize(BLMVM blmvm){
  CHKSOLVER(blmvm);
  blmvm->W=0;
  blmvm->DX=0;
  blmvm->GP=0;
  blmvm->G=0;
  blmvm->XL=0;
  blmvm->XU=0;
  blmvm->DXL=0;
  blmvm->DXU=0;
  blmvm->IHI=0;
  blmvm->M=0;
  blmvm->setup=0;
  blmvm->maxdx=1.0e200;
  blmvm->dxnorm=0;
  blmvm->alpha=0;
  return(0);
}

/* ---------------------------------------------------------- */
#undef __FUNCT__  
#define __FUNCT__ "BLMVMTakeDown"
static int BLMVMTakeDown(BLMVM blmvm){
  int info;
  CHKSOLVER(blmvm);
  if (blmvm->setup){
    info=BLMVMVecDestroy(&blmvm->W);CHKERRQ(info);
    info=BLMVMVecDestroy(&blmvm->DX);CHKERRQ(info);
    info=BLMVMVecDestroy(&blmvm->GP);CHKERRQ(info);
    info=BLMVMVecDestroy(&blmvm->G);CHKERRQ(info);
    info=BLMVMVecDestroy(&blmvm->XL);CHKERRQ(info);
    info=BLMVMVecDestroy(&blmvm->XU);CHKERRQ(info);
    info=BLMVMVecDestroy(&blmvm->DXL);CHKERRQ(info);
    info=BLMVMVecDestroy(&blmvm->DXU);CHKERRQ(info);
    info=LMVMMatDestroy(blmvm->M);CHKERRQ(info);
    info=BLMVMLineSearchTakeDown(blmvm->Line);CHKERRQ(info);
    info=BLMVMInitialize(blmvm);CHKERRQ(info);
  }
  return(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BLMVMGetStopReason"
static int BLMVMGetStopReason(BLMVM blmvm, int*flag){
  CHKSOLVER(blmvm);
  if (flag){*flag=blmvm->stopreason;}
  return(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BLMVMSetHistory"
static int BLMVMSetHistory(BLMVM blmvm, int lm){
  CHKSOLVER(blmvm);
  if (lm>=0){ blmvm->lm=lm;}
  return(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BLMVMGetHistory"
static int BLMVMGetHistory(BLMVM blmvm, int *lm){
  CHKSOLVER(blmvm);
  if (lm){*lm=blmvm->lm;}
  return(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BLMVMSetMaxStepNorm"
static int BLMVMSetMaxStepNorm(BLMVM blmvm, double stepnorm){
  CHKSOLVER(blmvm);
  if (stepnorm>0){ blmvm->maxdx=stepnorm;}
  return(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BLMVMGetStepNorm"
static int BLMVMGetStepNorm(BLMVM blmvm, double *stepnorm){
  CHKSOLVER(blmvm);
  if (stepnorm){ *stepnorm=blmvm->dxnorm;}
  return(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BLMVMGetStepSize"
static int BLMVMGetStepSize(BLMVM blmvm, double *stepsize){
  CHKSOLVER(blmvm);
  if (stepsize){ *stepsize=blmvm->alpha;}
  return(0);
}

#undef __FUNCT__  
#define __FUNCT__ "BLMVMGetSolutionVec"
static int BLMVMGetSolutionVec(BLMVM blmvm, BLMVMVec *X){
  CHKSOLVER(blmvm);
  if (X){*X=blmvm->X;}
  return(0);
}
#undef __FUNCT__  
#define __FUNCT__ "BLMVMGetProjectedGradientVec"
static int BLMVMGetProjectedGradientVec(BLMVM blmvm, BLMVMVec *PG){
  CHKSOLVER(blmvm);
  if (PG){*PG=blmvm->GP;}
  return(0);
}
#undef __FUNCT__  
#define __FUNCT__ "BLMVMGetGradientVec"
static int BLMVMGetGradientVec(BLMVM blmvm, BLMVMVec *G){
  CHKSOLVER(blmvm);
  if (G){*G=blmvm->G;}
  return(0);
}
#undef __FUNCT__  
#define __FUNCT__ "BLMVMGetStepVec"
static int BLMVMGetStepVec(BLMVM blmvm, BLMVMVec *DX){
  CHKSOLVER(blmvm);
  if (DX){*DX=blmvm->DX;}
  return(0);
}
#undef __FUNCT__  
#define __FUNCT__ "BLMVMGetBoundDualVec"
static int BLMVMGetBoundDualVec(BLMVM blmvm, BLMVMVec *DXL, BLMVMVec *DXU){
  CHKSOLVER(blmvm);
  if (DXL){*DXL=blmvm->DXL;}
  if (DXU){*DXU=blmvm->DXU;}
  return(0);
}
/* ----------------- End blmvm.c --------------- */
/* --------- Begin blmvmapplication.c ---------- */

#ifdef CHKAPP
#undef CHKAPP
#endif

#define CHKAPP(a) {if (!(a)||((a)->key!=AKEY)){ printf("BLMVMERROR: Invalid BLMVM_APP object\n  %s\n",__FUNCT__); return (1);}}

#ifdef AKEY
#undef AKEY
#endif
#define AKEY 2378


struct _P_BLMVM_APPLICATION{
  int(*fg)(void*, double*, double*, int, double*,double*); 
  int(*bounds)(void*, double* ,double* ,int);
  int (*converge)(void*, int,double,double*, int);
  void *appctx_fg;
  void *appctx_bounds;
  void *appctx_converge;
  int key;
};

#undef __FUNCT__  
#define __FUNCT__ "BLMVMAppCreate"
static int BLMVMAppCreate(BLMVM_APP *blmvmapp){
  BLMVM_APP blmvmapp2;
  BLMVMCALLOC(blmvmapp2,1,struct _P_BLMVM_APPLICATION);  
  blmvmapp2->key=AKEY;
  blmvmapp2->fg=0;
  blmvmapp2->bounds=0;
  blmvmapp2->converge=0;
  blmvmapp2->appctx_fg=0;
  blmvmapp2->appctx_bounds=0;
  blmvmapp2->appctx_converge=0;
  *blmvmapp=blmvmapp2;
  return(0);
}
#undef __FUNCT__  
#define __FUNCT__ "BLMVMAppDestroy"
static int BLMVMAppDestroy(BLMVM_APP *blmvmapp){
  CHKAPP(*blmvmapp);
  BLMVMFREE(*blmvmapp);
  *blmvmapp=0;
  return(0);
}

/* ---------------------------------------------------------- */

#undef __FUNCT__  
#define __FUNCT__ "BLMVMSetFunctionGradient"
static int BLMVMSetFunctionGradient(BLMVM_APP blmvmapp, int(*fg)(void*, double*, double*, int, double*, double*),void* ctx){
  CHKAPP(blmvmapp);
  blmvmapp->fg=fg;
  blmvmapp->appctx_fg=ctx;
  return (0);
}

#undef __FUNCT__  
#define __FUNCT__ "BLMVMSetBounds"
static int BLMVMSetBounds(BLMVM_APP blmvmapp, int(*bounds)(void*, double* ,double* ,int), void* ctx){
  CHKAPP(blmvmapp);
  blmvmapp->bounds=bounds;
  blmvmapp->appctx_bounds=ctx;
  return (0);
}


#undef __FUNCT__  
#define __FUNCT__ "BLMVMSetConverge"
static int BLMVMSetConverge(BLMVM_APP blmvmapp, int(*converge)(void*,int,double,double*, int), void* ctx){
  CHKAPP(blmvmapp);
  blmvmapp->converge=converge;
  blmvmapp->appctx_converge=ctx;
  return (0);
}

#undef __FUNCT__  
#define __FUNCT__ "BLMVMComputeFunctionGradient"
static int BLMVMComputeFunctionGradient(BLMVM_APP blmvmapp, BLMVMVec X,double *f,BLMVMVec G,BLMVMVec InvHessDiag){
  int info,n;
  double *x,*g, *ihd;
  CHKAPP(blmvmapp);
  info=BLMVMVecGetDoubles(X,&x,&n);CHKERRQ(info);
  info=BLMVMVecGetDoubles(G,&g,&n);CHKERRQ(info);
  info=BLMVMVecGetDoubles(InvHessDiag,&ihd,&n);CHKERRQ(info);
  if (blmvmapp->fg){
     info = (blmvmapp->fg)(blmvmapp->appctx_fg, x, g, n, f, ihd);CHKERRQ(info); 
  } else {
    printf("ERROR: NO FG Function pointer has been set\n");
    return 1;
  }
  return (0);
}

#undef __FUNCT__  
#define __FUNCT__ "BLMVMComputeBounds"
static int BLMVMComputeBounds(BLMVM_APP blmvmapp, BLMVMVec XL, BLMVMVec XU){
  int info,n;
  double *xl,*xu;
  CHKAPP(blmvmapp);
  info=BLMVMVecGetDoubles(XL,&xl,&n);CHKERRQ(info);
  info=BLMVMVecGetDoubles(XU,&xu,&n);CHKERRQ(info);
  if (blmvmapp->bounds){
      info=(blmvmapp->bounds)(blmvmapp->appctx_bounds, xl, xu, n);CHKERRQ(info);
  } else {
    printf("ERROR: NO Bounds function pointer has been set\n");
    return 1;
  }

  return (0);
}

#undef __FUNCT__  
#define __FUNCT__ "BLMVMMonitor"
static int BLMVMMonitor(BLMVM_APP blmvmapp, int iter, double step, double dxnorm, BLMVMVec PG){
  int n,info,finished;
  double *res;
  CHKAPP(blmvmapp);
  //printf("Going to GetDoubles\n");
  info=BLMVMVecGetDoubles(PG,&res,&n);CHKERRQ(info);
  //printf("Back from GetDoubles\n");
  finished=(blmvmapp->converge)(blmvmapp->appctx_converge,iter,step,res,n);
  return(finished);
}
/* --------- End blmvmapplication.c ---------- */
/* ------------ Begin blmvmmatrix.c ------------------------------- */
/* THE FOLLOWING CODE IMPLEMENTS APPROXIMATE HESSIAN INVERSE MATRIX */

#define min(a,b) ((a <= b)? (a) : (b))
#define max(a,b) ((a >= b)? (a) : (b))

#ifdef MKEY
#undef MKEY
#endif
#define MKEY 5284

#ifdef CHKMAT
#undef CHKMAT
#endif
#define CHKMAT(a) {if (!(a)||((a)->key!=MKEY)){ printf("BLMVMERROR: Invalid LMVMMat object\n  %s\n",__FUNCT__); return (1);}}

struct  _P_LMVMMat{

  int lm;
  int lmnow;
  int iter;
  int rejects;

  double eps;
  int key;
  BLMVMVec *S;
  BLMVMVec *Y;
  BLMVMVec Gprev;
  BLMVMVec Xprev;

  double y0normsquared;
  double *rho;
  double *beta;
};

int LMVMMatRefresh(LMVMMat M){
  CHKMAT(M);
  M->lmnow=0;
  M->iter=0;
  M->rejects=0;
  M->Gprev=0;
  M->Xprev=0;
  return (0);
}

int LMVMMatCreate(int nlm, BLMVMVec X, LMVMMat *MM){
  int i,info;
  LMVMMat M;
  BLMVMCALLOC(M,1,struct _P_LMVMMat);
  *MM=M;

  BLMVMCALLOC(M->S,(nlm+1),BLMVMVec);
  BLMVMCALLOC(M->Y,(nlm+1),BLMVMVec);

  BLMVMCALLOC(M->rho,(nlm+1),double);
  BLMVMCALLOC(M->beta,(nlm+1),double);
  
  M->lm=nlm;
  M->eps=2.2e-16;
  M->eps=2.2e-11;
  M->key=MKEY;

  for (i=0;i<nlm+1;i++){
    info=BLMVMVecDuplicate(X,&M->S[i]); CHKERRQ(info);
    info=BLMVMVecDuplicate(X,&M->Y[i]); CHKERRQ(info);
  }
  info= LMVMMatRefresh(M); CHKERRQ(info);
  return(0);
}

#undef __FUNCT__
#define __FUNCT__ "LMVMMatUpdate"
int LMVMMatUpdate(LMVMMat M, BLMVMVec  X, BLMVMVec G){
  int i,info;
  int lm=M->lm,lmnow=M->lmnow;
  double rhotemp,rhotol;
  double y0temp;
  double   *rho=M->rho;
  BLMVMVec *Y=M->Y;
  BLMVMVec *S=M->S;
  BLMVMVec Gprev=M->Gprev;
  BLMVMVec Xprev=M->Xprev;

  CHKMAT(M);
  if (M->Gprev==0 || M->Xprev==0){
    
    M->Gprev=M->Y[lm];
    M->Xprev=M->S[lm];
    M->rho[0]=1.0;
    M->y0normsquared = 1.0;
    M->iter=0;
    M->rejects=0;

  } else {
    
    M->iter++;
    info=BLMVMVecAYPX(-1.0,G,Gprev);CHKERRQ(info);
    info=BLMVMVecAYPX(-1.0,X,Xprev);CHKERRQ(info);
    info=BLMVMVecDot(Xprev,Gprev,&rhotemp);CHKERRQ(info);
    info=BLMVMVecDot(Gprev,Gprev,&y0temp);CHKERRQ(info);

    rhotol=M->eps*y0temp;
    if (rhotemp > rhotol){
		M->lmnow = min(lmnow+1, lm); 
//		M->lmnow = M->iter > 80 ? min(lmnow+1,lm) : min(lmnow+1, 40);
      for (i=lm-1;i>=0;i--){
        S[i+1]=S[i];
        Y[i+1]=Y[i];
        rho[i+1]=rho[i];
      }
      S[0]=M->Xprev;
      Y[0]=M->Gprev;
      rho[0]=1.0/rhotemp;
      M->y0normsquared=y0temp;
      M->Xprev=S[lm]; M->Gprev=Y[lm];

    } else { 
      M->rejects++;
    }
  }
  info=BLMVMVecCopy(X,M->Xprev);CHKERRQ(info);
  info=BLMVMVecCopy(G,M->Gprev);CHKERRQ(info);
  return (0);
}

#undef __FUNCT__
#define __FUNCT__ "LMVMMatSolve"
 int LMVMMatSolve(LMVMMat M, BLMVMVec G, BLMVMVec DX, BLMVMVec IHI){
   
  int      ll,info;
  double   sq, yq;
  double   *rho,*beta;
  BLMVMVec *Y,*S;
  double ans;
  double scale1=0, scale2=0;

  CHKMAT(M);
  Y=M->Y;S=M->S;
  rho=M->rho;beta=M->beta;

  if (M->lmnow<1){rho[0]=1.0; M->y0normsquared = 1.0;}

  info=BLMVMVecCopy(G,DX);CHKERRQ(info);  
  for (ll = 0; ll<M->lmnow; ll++){
    info=BLMVMVecDot(DX,S[ll],&sq);CHKERRQ(info);
    beta[ll] = sq*rho[ll];
    info=BLMVMVecAXPY(-beta[ll], Y[ll],DX);CHKERRQ(info);
  }

  // HERE
  info=BLMVMVecScale(1.0/(rho[0]*M->y0normsquared),DX);CHKERRQ(info);
  //info=BLMVMVecMultiply(-1.0, IHI, DX, 1.0/(rho[0]*M->y0normsquared), &scale1, &scale2);
  //printf("\n");
  //BLMVMVecView(IHI);
  //printf("\n");
  //BLMVMVecDot(IHI, DX, &ans);
  printf("   %0.20f %0.20f \n", (10.0/(rho[0]*M->y0normsquared)), scale2/scale1);

  for (ll=M->lmnow-1; ll>=0; ll--){
    info=BLMVMVecDot(DX,Y[ll],&yq);CHKERRQ(info);
    info=BLMVMVecAXPY(beta[ll]-yq*rho[ll],S[ll],DX);CHKERRQ(info);
  }
  return (0);
}

int LMVMMatX0(LMVMMat M, BLMVMVec X0){
  int i,info;
  CHKMAT(M);
  info=BLMVMVecCopy(M->Xprev,X0);CHKERRQ(info);
  for (i=0;i<M->lmnow;i++){
    info = BLMVMVecAXPY(-1.0,M->S[i],X0); CHKERRQ(info);
  }
  return (0);
}

int LMVMMatDestroy(LMVMMat M){
  int i,info;
  CHKMAT(M);
  for (i=0;i<M->lm+1;i++){
    info = BLMVMVecDestroy(&M->S[i]); CHKERRQ(info);
    info = BLMVMVecDestroy(&M->Y[i]); CHKERRQ(info);
  }
  M->key=0;
  BLMVMFREE(M->S);
  BLMVMFREE(M->Y);
  BLMVMFREE(M->rho);
  BLMVMFREE(M->beta);
  BLMVMFREE(M);
  return (0);
}
/* ------------ End blmvmmatrix.c ------------------------------- */



/* ---------- Begin blmvmlinesearch.c --------------------------- */


#ifdef LKEY
#undef LKEY
#endif
#define LKEY 9876

//#define min(a,b) ((a <= b)? (a) : (b))
//#define max(a,b) ((a >= b)? (a) : (b))

#define blmvmloginfo if(0)printf

struct  _P_BLMVMLineSearch{
  double backscale;
  int key;
};

#ifdef CHKLINE
#undef CHKLINE 
#endif
#define CHKLINE(a) {if (!(a)||((a)->key!=LKEY)){ printf("BLMVMERROR: %s\n  Invalid BLMVM_LINESEARCH object\n",__FUNCT__); return (1);}}

static int BLMVMStep_LineSearch(double*, double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

#undef __FUNCT__
#define __FUNCT__ "BLMVMApplyLineSearch"
static int BLMVMApplyLineSearch(BLMVM_LINESEARCH blmvmline, BLMVM_APP blmvmapp, BLMVMVec V[9],double *f, double*sstep){
  
  /*
   * Note: This routine expects positive step directions and tests for a valid descent direction of 
   *       DX accordingly. Because a negateve *step is passed in, DX and *step are 
   *       negated both in the beginning and at the end of this routine.
   */

  double    zero = 0.0, p5 = 0.5, p66 = 0.66, two =2.0, xtrapf = 4.0;
  double    finit, width, width1, dginit,fm, fxm, fym, dgm, dgxm, dgym;
  double    dgx, dgy, dg, fx, fy, stx, sty, dgtest, ftest1=0.0,ftest2=0.0;
  double    bstepmin1,bstepmin2,bstepmax;
  double    dg1,dg2,step;
  int       info,info2, i, stage1;
  int       infoc = 1;
  double stepmin  = 1.0e-30; /* lower bound for step */
  double stepmax  = 1.0e+20; /* upper bound for step */
  double rtol     = 1.0e-10; /* relative tolerance for an acceptable step */
  double ftol     = .05;     /* tolerance for sufficient decrease condition */
  double gtol     = .9;      /* tolerance for curvature condition */
  double nfev     = 0;       /* number of function evaluations */
  double maxfev   = 30;      /* maximum number of function evaluations */
  double bracket; 
  BLMVMVec XL=V[0],XU=V[1],DX=V[2],X=V[3],G=V[4],W=V[5],WA=V[6],WB=V[7],IHI=V[8];
  CHKLINE(blmvmline);
  
  /* The unconventional assumption is that the steplength is assumed to be negative */
  /*negate *step and DX to match standard for this routine (positive step) */
  step=-*sstep;
  printf("step: %f\n", step);

  /* Check input parameters for errors */
  if (step < zero) { info2 = -1; printf("ERROR: step < 0....Exiting \n"); return 0; } 
  else if (ftol < zero) { info2 = -2; printf("ERROR: ftol < 0....Exiting\n"); return 0;} 
  else if (rtol < zero) { info2 = -3; printf("ERROR: rtol < 0....Exiting\n"); return 0;} 
  else if (gtol < zero) { info2 = -4; printf("ERROR: gtol < 0....Exiting\n"); return 0;} 
  else if (stepmin < zero) { info2 = -5; printf("ERROR: stepmin < 0....Exiting\n"); return 0; } 
  else if (stepmax < stepmin) { info2 = -6; printf("ERROR: stepmax < stepmin....Exiting\n"); return 0;} 
  else if (maxfev < zero) { info2 = -7; printf("ERROR: makefev < 0....Exiting\n"); return 0;}

  /* Compute the smallest steplength that will make one nonbinding variable equal the bound */ 
  info = BLMVMVecProjectGradient(DX,XL,X,XU,WA,DX,WB); CHKERRQ(info);

  info = BLMVM_StepBound(X,XL,XU,DX,&bstepmin1,&bstepmin2,&bstepmax); CHKERRQ(info);
  if(stepmin > bstepmin1) { printf("BLMVMApplyLineSearch: min step too big (solution non-feasible)\n"); return 0; }
  stepmax = min(bstepmax,1.0e+15); 
  
  /* Check that search direction is a descent direction */
  info = BLMVMVecDot(G,DX,&dginit); CHKERRQ(info);
  dginit*=-1.0;
  if (dginit >= zero) { printf("Not a descent direction...Exiting"); info2 = 7; return 0;}

  /* Initialization */
  bracket=0; info2=0; stage1=1; finit=*f; dgtest=ftol*dginit; width=stepmax-stepmin; width1=width*two;
  info = BLMVMVecCopy(X,W); CHKERRQ(info);
  
  /* Variable dictionary:  
     stx, fx, dgx - the step, function, and derivative at the best step
     sty, fy, dgy - the step, function, and derivative at the other endpoint 
                    of the interval of uncertainty
     step, f, dg - the step, function, and derivative at the current step */
  
  stx = zero; fx  = finit; dgx = dginit; sty = zero; fy  = finit; dgy = dginit;
  nfev = 0;
  
  for (i=0; i< maxfev; i++) {
    /* Set min and max steps to correspond to the interval of uncertainty */
    if (bracket) { stepmin = min(stx,sty); stepmax = max(stx,sty); } 
    else { stepmin = stx; stepmax = step + xtrapf * (step - stx); }

    /* Force the step to be within the bounds */
    step = max(step,stepmin);
    step = min(step,stepmax);
    
    /* If an unusual termination is to occur, then let step be the lowest
       point obtained thus far */
    if (((bracket) && (step <= stepmin || step >= stepmax)) ||
        ((bracket) && (stepmax - stepmin <= rtol * stepmax)) ||
        (nfev >= maxfev - 1) || (infoc == 0)) 
        {
        step = stx;
        }
    
	printf("step here is %f\n", step);
//	printf("X is ");
//	BLMVMVecView(X);
//	printf("DX is ");
//	BLMVMVecView(DX);

    info = BLMVMVecWAXPY(-step,DX,W,X); CHKERRQ(info); 	/* X = W + step*DX */
    
    /* This if statement added by S. Benson.   --  Project the solution, if necessary, to keep feasible */
    if (step >= bstepmin1){
      info=BLMVMVecApplyBounds(X,XL,XU); CHKERRQ(info);
    }
//	if (step < 1E-4) {
//		break;
//	}
    
    info = BLMVMComputeFunctionGradient(blmvmapp,X,f,G,IHI); CHKERRQ(info);
//	printf("Computed gradient: ");
//	BLMVMVecView(G);
    nfev++;
    
    /*    info = G->Dot(S,&dg);CHKERRQ(info);	        / * dg = G^T S */
     /* Next 3 lines  added by S. Benson.   to compute an estimate of decrease in the actual direction and curvature */
    info = BLMVMVecDot(G,W,&dg1); CHKERRQ(info);	        /* dg = G^T S */
    info = BLMVMVecDot(G,X,&dg2); CHKERRQ(info);	        /* dg = G^T S */
    dg=dg2-dg1;
    ftest1 = finit + (step) * dgtest;

    /* This test added by S. Benson.   --  Be satisfied with Armijo condition if the Projection X->Median(XL,X,XU) changed X. */
    ftest2 = finit + (step) * dgtest * ftol;
    
    /* Convergence testing */
    if (((bracket) && (step <= stepmin||step >= stepmax)) || (!infoc)) 
	{ 
		printf("BLMVMApply_BoundLineSearch:Possible Rounding errors.\n"); 
	info2 = 6; 
	break; 
	}
    if ((step == stepmax) && (*f <= ftest1) && (dg <= dgtest)) 
    { printf("BLMVMApply_BoundLineSearch:Step is at the upper bound, stepmax (%g)\n",stepmax); info2 = 5; break;}
    if ((step == stepmin) && (*f >= ftest1) && (dg >= dgtest)) 
    { printf("BLMVMApply_BoundLineSearch:Step is at the lower bound, stepmin (%g)\n",stepmin); info2 = 4; break;}
    if (nfev >= maxfev) 
    { printf("BLMVMApply_BoundLineSearch:Number of line search function evals  > maximum \n"); info2 = 3; break; }
    if ((bracket) && (stepmax - stepmin <= rtol*stepmax)) 
    { printf("BLMVMApply_BoundLineSearch:Relative width of interval of uncertainty is at most rtol (%g)\n",rtol); info2 = 2; break;}
    if ((*f <= ftest1) && (fabs(dg) <= gtol*(-dginit))) 
    { /*success*/info2 = 0; 
//	printf("succeeded in inner loop; not doing line search.  %f %f, %f, %f", (*f), (ftest1), (fabs(dg)), (gtol*(-dginit)));
	break; }

    /* This test added by S. Benson.   --  Be satisfied with Armijo condition if the Projection X->Median(XL,X,XU) changed X. */
    if ( (*f <= ftest2) && ((step) < bstepmin2) ) 
    {blmvmloginfo("BLMVMApply_BoundLineSearch:Line search success: Sufficient decrease new bounds apply\n"); info2 = 0; break;}

    /* In the first stage, we seek a step for which the modified function
        has a nonpositive value and nonnegative derivative */
    if ((stage1) && (*f <= ftest1) && (dg >= dginit * min(ftol, gtol))) stage1 = 0;

    /* A modified function is used to predict the step only if we
       have not obtained a step for which the modified function has a 
       nonpositive function value and nonnegative derivative, and if a
       lower function value has been obtained but the decrease is not
       sufficient */

    if ((stage1) && (*f <= fx) && (*f > ftest1)) {
      fm   = *f - step * dgtest;	/* Define modified function */
      fxm  = fx - stx * dgtest;	        /* and derivatives */
      fym  = fy - sty * dgtest;
      dgm  = dg - dgtest;
      dgxm = dgx - dgtest;
      dgym = dgy - dgtest;

      /* Update the interval of uncertainty and compute the new step */
      info = BLMVMStep_LineSearch(&stepmin,&stepmax,&bracket,&stx,&fxm,&dgxm,&sty,&fym,&dgym,&step,&fm,&dgm); CHKERRQ(info);

      fx  = fxm + stx * dgtest;	/* Reset the function and */
      fy  = fym + sty * dgtest;	/* gradient values */
      dgx = dgxm + dgtest; 
      dgy = dgym + dgtest; 
    } else {
      /* Update the interval of uncertainty and compute the new step */
      info = BLMVMStep_LineSearch(&stepmin,&stepmax,&bracket,&stx,&fx,&dgx,&sty,&fy,&dgy,&step,f,&dg); CHKERRQ(info);
    }

   /* Force a sufficient decrease in the interval of uncertainty */
   if (bracket) {
     if (fabs(sty - stx) >= p66 * width1) step = stx + p5*(sty - stx);
       width1 = width;
       width = fabs(sty - stx);
     }
   }

  /* Finish computations */
  /* printf("BLMVM_BoundLineSearch: step = %4.4e\n",step); */
  /*
  info = G->Norm2(dginitt);CHKERRQ(info);
  */

  /*reset to original standard*/ 
  *sstep=(-1)*(step);

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "BLMVMStep_LineSearch"
static int BLMVMStep_LineSearch(double *stepmin, double *stepmax,double *bracket,double *stx,double *fx,double *dx,
				double *sty,double *fy,double *dy,double *stp,
				double *fp,double *dp)
{
  double            gamma1, p, q, r, s, sgnd, stpc, stpf, stpq, theta;
  double            zero = 0.0;
  int               bound;
  int infoc;
  
  /* Check the input parameters for errors */
  infoc = 0;
  if (*bracket && (*stp <= min(*stx,*sty) || (*stp >= max(*stx,*sty))))
    {printf("BLMVMStep_LineSearch: bad step in bracket\n");}
  if (*dx * (*stp-*stx) >= zero) printf("BLMVMStep_LineSearch: dx * (stp-stx) >= 0\n");
  if (*stepmax < *stepmin) printf("BLMVMStep_LineSearch: stepmax > stepmin\n");

  /* Determine if the derivatives have opposite sign */
  sgnd = *dp * (*dx/fabs(*dx));

  /* Case 1: a higher function value.
     the minimum is bracketed. if the cubic step is closer
     to stx than the quadratic step, the cubic step is taken,
     else the average of the cubic and quadratic steps is taken.
   */
  if (*fp > *fx) {
    printf("case 1\n");
    infoc = 1;
    bound = 1;
    theta = 3 * (*fx - *fp) / (*stp - *stx) + *dx + *dp;
    s = max(fabs(theta),fabs(*dx));
    s = max(s,fabs(*dp));
    gamma1 = s*sqrt(pow(theta/s,2) - (*dx/s)*(*dp/s));
    if (*stp < *stx) gamma1 = -gamma1;
    /* Can p be 0?  Check */
    p = (gamma1 - *dx) + theta;
    q = ((gamma1 - *dx) + gamma1) + *dp;
    r = p/q;
    stpc = *stx + r*(*stp - *stx);
    stpq = *stx + ((*dx/((*fx-*fp)/(*stp-*stx)+*dx))*0.5) * (*stp - *stx);
    if (fabs(stpc-*stx) < fabs(stpq-*stx)) { stpf = stpc;} 
    else { stpf = stpc + 0.5*(stpq - stpc); }
    *bracket = 1;
  }
  /* 
     Case 2: A lower function value and derivatives of
     opposite sign. the minimum is bracketed. if the cubic
     step is closer to stx than the quadratic (secant) step,
     the cubic step is taken, else the quadratic step is taken.
  */
  else if (sgnd < zero) {
    printf("case 2\n");
    infoc = 2;
    bound = 0;
    theta = 3*(*fx - *fp)/(*stp - *stx) + *dx + *dp;
    s = max(fabs(theta),fabs(*dx));
    s = max(s,fabs(*dp));
    gamma1 = s*sqrt(pow(theta/s,2) - (*dx/s)*(*dp/s));
    if (*stp > *stx) gamma1 = -gamma1;
    p = (gamma1 - *dp) + theta;
    q = ((gamma1 - *dp) + gamma1) + *dx;
    r = p/q;
    stpc = *stp + r*(*stx - *stp);
    stpq = *stp + (*dp/(*dp-*dx))*(*stx - *stp);
    if (fabs(stpc-*stp) > fabs(stpq-*stp)) stpf = stpc;
    else stpf = stpq;
    *bracket = 1;
  }

/*   Case 3: A lower function value, derivatives of the
     same sign, and the magnitude of the derivative decreases.
     the cubic step is only used if the cubic tends to infinity
     in the direction of the step or if the minimum of the cubic
     is beyond stp. otherwise the cubic step is defined to be
     either stepmin or stepmax. the quadratic (secant) step is also
     computed and if the minimum is bracketed then the the step
     closest to stx is taken, else the step farthest away is taken.
 */

  else if (fabs(*dp) < fabs(*dx)) {
    printf("case 3\n");

    infoc = 3;
    bound = 1;
    theta = 3*(*fx - *fp)/(*stp - *stx) + *dx + *dp;
    s = max(fabs(theta),fabs(*dx));
    s = max(s,fabs(*dp));

    /* The case gamma1 = 0 only arises if the cubic does not tend
       to infinity in the direction of the step. */
    gamma1 = s*sqrt(max(zero,pow(theta/s,2) - (*dx/s)*(*dp/s)));
    if (*stp > *stx) gamma1 = -gamma1;
    p = (gamma1 - *dp) + theta;
    q = (gamma1 + (*dx - *dp)) + gamma1;
    r = p/q;
    if (r < zero && gamma1 != zero) stpc = *stp + r*(*stx - *stp);
    else if (*stp > *stx)        stpc = *stepmax;
    else stpc = *stepmin;
    stpq = *stp + (*dp/(*dp-*dx)) * (*stx - *stp);
    if (*bracket) {
      if (fabs(*stp-stpc) < fabs(*stp-stpq)) stpf = stpc;
      else stpf = stpq;
    }
    else {
      if (fabs(*stp-stpc) > fabs(*stp-stpq)) stpf = stpc;
      else stpf = stpq;
    }
  }

/*   Case 4: A lower function value, derivatives of the
     same sign, and the magnitude of the derivative does
     not decrease. if the minimum is not bracketed, the step
     is either stpmin or stpmax, else the cubic step is taken.
 */
  else {
    printf("case 4\n");

    infoc = 4;
    bound = 0;
    if (*bracket) {
      theta = 3*(*fp - *fy)/(*sty - *stp) + *dy + *dp;
      s = max(fabs(theta),fabs(*dy));
      s = max(s,fabs(*dp));
      gamma1 = s*sqrt(pow(theta/s,2) - (*dy/s)*(*dp/s));
      if (*stp > *sty) gamma1 = -gamma1;
      p = (gamma1 - *dp) + theta;
      q = ((gamma1 - *dp) + gamma1) + *dy;
      r = p/q;
      stpc = *stp + r*(*sty - *stp);
      stpf = stpc;
    } 
    else if (*stp > *stx) { stpf = *stepmax; } 
    else { stpf = *stepmin; }
  }

  /* Update the interval of uncertainty.  This update does not
     depend on the new step or the case analysis above. */

  if (*fp > *fx) { *sty = *stp; *fy = *fp; *dy = *dp;} 
  else {
    if (sgnd < zero) { *sty = *stx; *fy = *fx; *dy = *dx;}
    *stx = *stp;
    *fx = *fp;
    *dx = *dp;
  }

  /* Compute the new step and safeguard it */
  stpf = min(*stepmax,stpf);
  stpf = max(*stepmin,stpf);
  *stp = stpf;
  if (*bracket && bound) {
	  if (*sty > *stx) *stp = min(*stx+0.66*(*sty-*stx),*stp);
	  else             *stp = max(*stx+0.66*(*sty-*stx),*stp);
  }
  return 0;

}

#undef __FUNCT__
#define __FUNCT__ "BLMVMLineSearchCreate"
static int BLMVMLineSearchCreate(BLMVM_LINESEARCH *pblmvmline){
  BLMVMCALLOC(*pblmvmline,1,struct _P_BLMVMLineSearch);
  (*pblmvmline)->key=LKEY;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "BLMVMLineSearchSetUp"
static int BLMVMLineSearchSetUp(BLMVM_LINESEARCH blmvmline, BLMVMVec X){
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "BLMVMLineSearchTakeDown"
static int BLMVMLineSearchTakeDown(BLMVM_LINESEARCH blmvmline){
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "BLMVMLineSearchDestroy"
static int BLMVMLineSearchDestroy(BLMVM_LINESEARCH *blmvmline){
  BLMVMFREE(*blmvmline);
  return 0;
}

/* ------------- End blmvmlinesearch.c ----------- */
 
/* ----------- Begin blmvmvec.c ---------- */

//#define min(a,b) ((a <= b)? (a) : (b))
//#define max(a,b) ((a >= b)? (a) : (b))


#ifdef VKEY
#undef VKEY
#endif

#define VKEY 2345

struct  _P_BLMVMVec{
  int    dim;
  double *val;
  int key;
};


#ifdef CHKVEC
#undef CHKVEC
#endif
#define CHKVEC(a,b) {if (!(a)||((a)->key!=VKEY)){ printf("BLMVMERROR: %s\n  Invalid Argument: %d, Expecting BLMVMVec object\n",__FUNCT__,b); return (1);}}


/* THE FOLLOWING CODE IMPLEMENTS VECTOR OPERATIONS */
/* ---------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "BLMVMVecGetDoubles"
static int BLMVMVecGetDoubles(BLMVMVec V, double **dd, int *nn){
  CHKVEC(V,1);
  *dd=V->val;
  *nn=V->dim;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "BLMVM_StepBound"
/*compute smallest step which will bring at least one variable "out of bounds" (bboundmin)*/
/*compute smallest step which will bring all variables "out of bounds" (bboundmax)*/
static int BLMVM_StepBound(BLMVMVec X, BLMVMVec XL, BLMVMVec XU,BLMVMVec DX, double *boundmin,double *wolfemin,double *boundmax){
  int i,n;
  double *x, *xl, *xu, *dx;
  double tt,t1=1.0e+20, t2=1.0e+20, t3=0;
  double tt1,tt2,tt3;
      
  CHKVEC(XL,1); CHKVEC(XU,2); CHKVEC(X,3); CHKVEC(DX,4);
  
  n=X->dim;
  x=X->val;  xl=XL->val;  xu=XU->val;  dx=DX->val;

  for (i=0;i<n;i++){
    if (-dx[i]>0){
      tt=(xu[i]-x[i])/-dx[i];
      t1=min(t1,tt);
      if (t1>0){
	t2=min(t2,tt);
      }
      t3=max(t3,tt);
    } else if (-dx[i]<0) {
      tt=(xl[i]-x[i])/-dx[i];
      t1=min(t1,tt);
      if (t1>0){
	t2=min(t2,tt);
      }
      t3=max(t3,tt);
    }
  }
  tt1=t1; tt2=t2; tt3=t3;
  *boundmin=tt1;
  *wolfemin=tt2;
  *boundmax=tt3;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "BLMVMVecCreate"
static int BLMVMVecCreate(int n ,BLMVMVec *VV){
  BLMVMVec V;
  BLMVMCALLOC(V,1,struct _P_BLMVMVec);
  V->dim=n;
  V->key=VKEY;
  *VV=V;
  if (n>0){
    BLMVMCALLOC(V->val,n,double);
    if (V->val==0) return 1;
  } else {
    V->val=0;
  }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "BLMVMVecCreateWArray"
static int BLMVMVecCreateWArray(BLMVMVec *VV, double* vv, int n){
  BLMVMVec V;
  BLMVMCALLOC(V,1,struct _P_BLMVMVec);
  V->dim=n;
  V->key=VKEY;
  *VV=V;
  if (n>0){
    V->val=vv;
  } else {
    V->val=0;
  }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "BLMVMVecDestroy"
static int BLMVMVecDestroy(BLMVMVec *V){
  CHKVEC(*V,1);
  
  if ((*V)->val){ BLMVMFREE((*V)->val); }
  BLMVMFREE(*V);
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "BLMVMVecCopy"
static int BLMVMVecCopy( BLMVMVec v1,  BLMVMVec v2){
  int n;
  double *val1, *val2;
  CHKVEC(v1,1); CHKVEC(v2,2);  
  n=v1->dim;  val1=v1->val; val2=v2->val;
  memcpy(val2,val1,n*sizeof(double));
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "BLMVMVecScale"
static int BLMVMVecScale(double alpha, BLMVMVec x){
  int i,n;
  double *xx;
  CHKVEC(x,2);
  n=x->dim;  xx=x->val;
  for (i=0; i<n; ++i){ xx[i]*= alpha;}
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "BLMVMVecAXPY"
static int BLMVMVecAXPY(double alpha,  BLMVMVec x,  BLMVMVec y){
  int i,n;
  double *yy, *xx;
  CHKVEC(x,2); CHKVEC(y,3);
  n=x->dim;  yy=y->val;  xx=x->val;
  for (i=0; i<n; ++i){ yy[i] += (alpha)*xx[i];}
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "BLMVMVecMultiply"
static int BLMVMVecMultiply(double alpha,  BLMVMVec x,  BLMVMVec y, double limit, double* measure1, double* measure2){
  int i,n;
  double *yy, *xx;
  CHKVEC(x,2); CHKVEC(y,3);
  n=x->dim;  yy=y->val;  xx=x->val;

  *measure1 = 0;
  *measure2 = 0;
  for (i=0; i<n; ++i){ 
	  double mult = alpha*xx[i];
	  if (mult>0 && mult > limit) mult = limit;
	  if (mult<0 && -mult>-limit) mult = -limit;
	  *measure1 += fabs(yy[i]);
	  yy[i] *= mult;
	  *measure2 += fabs(yy[i]);
	  //printf("%f %f %f %f\n", *measure1, *measure2, yy[i], mult);
  }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "BLMVMVecWAXPY"
static int BLMVMVecWAXPY(double alpha, BLMVMVec x, BLMVMVec y, BLMVMVec w){
  int i,n;
  double *yy, *xx, *ww;
  CHKVEC(x,2); CHKVEC(y,3); CHKVEC(w,4);  
  n=x->dim, yy=y->val, xx=x->val, ww=w->val;
  for (i=0; i<n; i++){ ww[i] = yy[i]+(alpha)*xx[i]; }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "BLMVMVecView"
static int BLMVMVecView(BLMVMVec V){
  int i,n;
  double *vv;
  CHKVEC(V,1);
  n=V->dim;  vv=V->val;
  for (i=0; i<n; ++i){printf("%4.4e ",vv[i]);} printf("\n");  
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "BLMVMVecAYPX"
static int BLMVMVecAYPX(double alpha,  BLMVMVec x,  BLMVMVec y){
  int i,n;
  double *yy, *xx;
  CHKVEC(x,2); CHKVEC(y,3);
  yy=y->val; xx=x->val;  n=x->dim;
  for (i=0; i<n; ++i){ yy[i]=xx[i]+(alpha)*yy[i]; }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "BLMVMVecDuplicate"
static int BLMVMVecDuplicate(BLMVMVec V1,BLMVMVec *V2){
  int info,n;
  CHKVEC(V1,1);  
  n=V1->dim;
  info = BLMVMVecCreate(n ,V2);
  return info;
}

#undef __FUNCT__
#define __FUNCT__ "BLMVMVecDot"
static int BLMVMVecDot(BLMVMVec V1, BLMVMVec V2, double *ans){
  int i,m;
  double *v1, *v2;
  CHKVEC(V1,1); CHKVEC(V2,2);  
  m=V1->dim;  v1=V1->val;  v2=V2->val;
  *ans=0.0;
  for (i=0; i<m; ++i){ *ans += v1[i]*v2[i]; }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "BLMVMVec2Norm"
static int BLMVMVec2Norm(BLMVMVec V, double *ans){
  int i,m;
  double *v;
  CHKVEC(V,1);
  m=V->dim;  v=V->val;
  *ans=0.0;
  for (i=0; i<m; ++i){ *ans += v[i]*v[i]; }
  *ans=sqrt(*ans);
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "BLMVMVecApplyBounds"
static int BLMVMVecApplyBounds( BLMVMVec V, BLMVMVec VMin, BLMVMVec VMax){
  int i,n=V->dim;
  double *v1,*v2,*v3;
  CHKVEC(V,1); CHKVEC(VMin,2); CHKVEC(VMax,3);
  n=V->dim; v1=V->val; v2=VMin->val; v3=VMax->val;
  for (i=0; i<n; ++i){ 
    if (v2[i]>v3[i]){ return 1;}
    v1[i]=max(v1[i],v2[i]); v1[i]=min(v3[i],v1[i]); 
  }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "BLMVMVecProjectGradient"
static int BLMVMVecProjectGradient( BLMVMVec G, BLMVMVec XL, BLMVMVec X, BLMVMVec XU, BLMVMVec DXL, BLMVMVec GP, BLMVMVec DXU){  
  int i,n;
  double *x,*xl,*xu,*g,*gp,*dxl,*dxu;
  CHKVEC(G,1);CHKVEC(XL,2); CHKVEC(X,3); CHKVEC(XU,4); CHKVEC(DXL,5); CHKVEC(GP,4); CHKVEC(DXU,7);
  x=X->val;  xl=XL->val;  xu=XU->val;  g=G->val; gp=GP->val; dxl=DXL->val; dxu=DXU->val;
  n=X->dim;
  
  for (i=0; i<n; i++){ 
    if (g[i]>0 && x[i]<=xl[i]){
      gp[i] = 0; dxl[i]=g[i]; dxu[i]=0;
    } else if (g[i]<0 && x[i]>=xu[i]){
      gp[i] = 0; dxl[i]=0; dxu[i]=g[i];
    } else {
      gp[i] = g[i]; dxl[i]=0; dxu[i]=0;
    }
  }
  return (0);
}

#undef __FUNCT__
#define __FUNCT__ "BLMVMVecCompatible"
static int BLMVMVecCompatible(BLMVMVec V1,BLMVMVec V2, int *flag){
  CHKVEC(V1,1); CHKVEC(V2,2);
  *flag=0;
  if (V1->dim==V2->dim){
    *flag=1;
  }
  return (0);
}

static int BLMVMVecDim(BLMVMVec v) {
	return v->dim;
}
/* ----------- End blmvmvec.c ---------- */
