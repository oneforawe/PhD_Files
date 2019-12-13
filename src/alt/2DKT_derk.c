/*  File Comments  *//*=========================================================================*//* FILENAME: vlt_derk.c *//* AUTHOR: Andrew Forrester <aforrester@ucla.edu>  *//* DESCRIPTION:   This C program employs the 4th-order Runge Kutta method in conjunction with the other vlt (vortex loop theory) programs (vlt_ThermStates.c, etc.) to solve for thermodynamic quantities as functions of temperature, length scale, and sometimes pressure.   In the filename, "vlt" stands for Vortex Loop Theory (of the superfluid phase transition) and "derk" stands for the Differential Equation Runge Kutta method.   See below the code for explanation of constants and variables.   Input:   Output:*//* REF FILE: vlt_ThermStates.c, vlt_K0cFind.c, vlt_HeatCap.c, ring1.c, 2DTK.c, etc. *//*  Function Preparation  *//*=========================================================================*//* Standard routine header files */#include <stdio.h>/* Constants and labels */#define IN    // input label#define OUT   // output label#define INOUT // input/output label/* Program parameters/inputs */#define N     2  // If using ring1.c N=3, If using 2DKT.c or K0cFind1.c N=2/* Data type definitions */typedef struct {	double l;	double tempv;	double thrmv[N];} STATE;typedef struct {	double arr[N];} RETARRAY;/*  Function definition  *//*=========================================================================*///   derk(func(),calculated,Dl,n)//   func(l,z)STATE derk(RETARRAY (*func)(double l, double z[]), STATE calculated, double Dl, unsigned n){	double dl;	int i,j,k;	RETARRAY DthrmvDl[4];	STATE test[3];	/* Calculate the derivatives at the input calculated state and 3n test states */	DthrmvDl[0] = func(IN calculated.l, IN calculated.thrmv);	for(i=0; i<3; i++){		dl = Dl*((i+2)/2)/2;  // See "Note" below.		for(j=0; j<n; j++){			test[i].thrmv[j] = calculated.thrmv[j] + DthrmvDl[i].arr[j]*dl;			test[i].l = calculated.l + Dl*((i+1)/2)/2;  // See "Note" below.		}		DthrmvDl[i+1] = func(IN test[i].l, IN test[i].thrmv);	}	// Note: This algorithm, using '((...)/2)/2', creates a desired rounding "error" that produces the sequence 0.0, 0.5, 0.5, 1.0, instead of 0.0, 0.5, 0.75, 1.0.	/* Update the calculated state, using the derivatives calculated above */	for(j=0; j<n; j++) calculated.thrmv[j] += Dl*(DthrmvDl[0].arr[j]+2*(DthrmvDl[1].arr[j]+DthrmvDl[2].arr[j])+DthrmvDl[3].arr[j])/6;	calculated.l += Dl;	/* Return the new calculated state */	return calculated;}/*  Program Notes  *//*=========================================================================*//*See further documentation (UCLAresearchNotes.pdf, UCLAresearchPrograms.pdf, vlt_ThermStates.c, ring1.c) for elaboration.Dl is the increment in l.n is the number of thermodynamic variables (n=3 given thermv[i], i=1,2,3).DthrmvDl[i].arr[0] = dK/dlDthrmvDl[i].arr[1] = dy/dlDthrmvDl[i].arr[2] = de/dl*/