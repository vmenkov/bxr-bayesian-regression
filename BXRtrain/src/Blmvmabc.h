/* --- Begin blmvmabc.h ---------- */
#ifndef BLMVMABC_H
#define BLMVMABC_H
extern int BLMVMFunctionGradient(void*,double[],double[],int,double*,double*);
extern int BLMVMConverge(void*,int,double,double[],int);
extern int BLMVMSolveIt(void*,double[],double[],double[],int,int);
extern int BLMVMSolve(int(*)(void*,double[],double[],int,double*),
		      int (*)(void*,int,double,double[],int),
		      void*,double[],double[],double[],int,int);
#endif
/* --- End blmvmabc.h ---------- */
