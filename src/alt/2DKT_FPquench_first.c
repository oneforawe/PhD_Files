/*  File Comments  */
/*=========================================================================*/

/* FILENAME: 2DKT_FPquench_first.c */
/* DESCRIPTION: 
 */
/* EXT FILES: 2DKT_derk.c, 2DKT_derk1.c */
/* AUTHOR: Andrew Forrester <aforrester@ucla.edu> */



/*  Function Preparation  */
/*=========================================================================*/

/* Standard routine header files */
#include <stdio.h>
#include <math.h>

/* Constants and labels */
#define IN      // input label
#define OUT     // output label
#define INOUT   // input/output label
#define PI      3.14159265358979323846
#define PISQ    9.86960440108935799230
//              9.86960440108935861883 ?
#define PICU    31.00627668029981620634
#define B       4.0*PICU

/* Program parameters/inputs */
#define N       2
// #define Kthresh 5.0
#define a0      1
#define a04     1  //  a0 to the fourth power
#define K0c     0.747853
#define dl      0.0001
#define lmax    10
// lpts = lmax/dl = 10/0.0001 = 100,000  (make sure this is divisible by datamax)
// the number of "l" points will really be lpts+1
#define tumax   5      //  (later, 10000) max unitless time (see if(tu==...) below)
#define iTfrac  0.95   //  initial temperature fraction Ti/Tkt
#define qTfrac  0.10   //  quench temperature fraction Tq/Tkt
#define datamax 100    //  max number of data points per "snap-shot" in time

/* Data type definitions */
typedef struct {
	double l;
	double tempv;
	double thrmv[N];
} STATE;

typedef struct {
	double arr[N];
} RETARRAY;



/*  Function prototypes  */
/*=========================================================================*/
//   derk(func(),calculated,Dl,n)
//   derk1(func(),calculated,Dl,n)
//   VPairRecRel(l,z)
//   FPRecRel(l,z)

extern STATE derk(RETARRAY (*func)(double l, double z[]), STATE calculated, double Dl, unsigned n);

extern double derk1(RETARRAY (*func)(double l, double z[]), STATE calculated, double Dl, unsigned n);

RETARRAY VPairRecRel(double l, double z[]);

RETARRAY FPRecRel(double l, double z[]);


/*  Function definitions  */
/*=========================================================================*/

main(){
	int tu,i,lpts;
	double K0,y0;
	lpts = lmax/dl;
	double Geq[lpts+1], dGdt[lpts+1];
	STATE estimated[lpts+1];
	FILE *outfile;
	char filename[100];

	/* Print headings for data to screen and output file */
	sprintf(filename, "2DKT_FPquench_lmax%i_K0c%7.6f.dat", lmax, K0c);
	outfile = fopen(filename,"w");  //  E.g., "2DKT_FPquench.out" or "2DKT_FPquench_lmax100_K0c0.747853.dat"
	fprintf(outfile,"%s\t%s\t%s\t%s\t%s\t%s\n",
		"T/Tc","1-T/Tc","l","K","G","K/K0");
	printf(         "%s\t%s\t%s\t%s\t%s\t%s\n",
		"T/Tc","1-T/Tc","l","K","G","K/K0");


    /*** Initial (tu=0), Equilibrium Calculations ***/


	/* Set initial temperature (Tkt) */
	estimated[0].tempv = 1.0-iTfrac;

	/* Initialize length scale value */
	estimated[0].l = 0.0;

	/* Initialize temperature-dependent values at initial length scale */
	K0   = K0c/(1.0-estimated[0].tempv);
	estimated[0].thrmv[0] = K0;               //  K
	estimated[0].thrmv[1] = exp(-PISQ*K0/2);  //  y
	Geq[0] = exp(-4*estimated[0].l)*estimated[0].thrmv[1]*estimated[0].thrmv[1]/a04;  // not = 0.0;

	/* Print initial data to screen and output file */
	fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\n",
		1-estimated[0].tempv, estimated[0].tempv, estimated[0].l, estimated[0].thrmv[0], Geq[0], estimated[0].thrmv[0]/K0);  // writes to a temporary file (not yet to outfile, without fflush)
	printf(         "%g\t%g\t%g\t%g\t%g\t%g\n",
		1-estimated[0].tempv, estimated[0].tempv, estimated[0].l, estimated[0].thrmv[0], Geq[0], estimated[0].thrmv[0]/K0);

	/* Calculate state-variables at progressively larger length scales, up to coherence length, via Runge-Kutta method (derk) */
	for(i=1; i<=lpts; i++){
		// if(estimated[i].thrmv[0]>Kthresh) break;

		/* Calculate next estimated state, incrementing length scale */
		estimated[i] = derk(IN VPairRecRel, IN estimated[i-1], IN dl, IN N);
		estimated[i].tempv = 1.0-iTfrac;
		Geq[i] = exp(-4*estimated[i].l)*estimated[i].thrmv[1]*estimated[i].thrmv[1]/a04;

		/* Print at most [datamax] (e.g., 1000) data points to screen and output file */
		if(i%(lpts/datamax)==0){
			fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\n",
				1-estimated[i].tempv, estimated[i].tempv, estimated[i].l, estimated[i].thrmv[0], Geq[i], estimated[i].thrmv[0]/K0);  // writes to a temporary file (not yet to outfile, without fflush)
			printf(         "%g\t%g\t%g\t%g\t%g\t%g\n",
				1-estimated[i].tempv, estimated[i].tempv, estimated[i].l, estimated[i].thrmv[0], Geq[i], estimated[i].thrmv[0]/K0);
		}
	}

	/* For convenience, set y as G */
	for(i=0; i<=lpts; i++){
		estimated[i].thrmv[1] = Geq[i];
	}


    /*** After Temperature Drop/Quench ***/


	/* Advance in time */
	for(tu=1; tu<=tumax; tu++){

		/* Set initial bath temperature (Tq) */
		estimated[0].tempv = 1.0-qTfrac;

		/* Initialize length scale value */
		// estimated[0].l = 0.0; // (already done)

		/* Initialize temperature-dependent values at initial length scale */
		K0   = K0c/(1.0-estimated[0].tempv);
		estimated[0].thrmv[0] = K0;  //  K
		y0 = exp(-PISQ*K0/2);        //  y
		estimated[0].thrmv[1] = exp(-4*estimated[0].l)*y0*y0/a04;  //  G  (not = 0.0;)

		/* At certain times, print out the results, starting, here, at the smallest length scale */
		if(tu<5){   //  tu==1||tu==10||tu==100||tu==1000||tu==10000){  //  (log10(tu)%1==0){
			/* Print initial data to screen and output file */
			fprintf(outfile,"\n%s%i\n\n", "tu = ",tu);
			fprintf(outfile,"%s\t%s\t%s\t%s\t%s\t%s\n",
				"T/Tc","1-T/Tc","l","K","G","K/K0");
			fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\n",
				1-estimated[0].tempv, estimated[0].tempv, estimated[0].l, estimated[0].thrmv[0], estimated[0].thrmv[1], estimated[0].thrmv[0]/K0);  // writes to a temporary file (not yet to outfile, without fflush)
			printf(         "\n%s%i\n\n%g\t%g\t%g\t%g\t%g\t%g\n", "tu = ",tu,
				1-estimated[0].tempv, estimated[0].tempv, estimated[0].l, estimated[0].thrmv[0], estimated[0].thrmv[1], estimated[0].thrmv[0]/K0);
		}

		/* Calculate dGdt, for all length scales */
		// dGdt[0] = 0; // since G[0] = 0, always
		for(i=1; i<lpts; i++){
			// if(estimated[i].thrmv[0]>Kthresh) break;
			//  dGdt[i] = exp(-2l)*( PI*(K[i+1]*G[i+1]-K[i-1]*G[i-1])/dl + (G[i+1]-2*G[i]+G[i-1]])/(dl*dl) );
			dGdt[i] = exp(-2*estimated[i].l)*( PI*(estimated[i+1].thrmv[0]*estimated[i+1].thrmv[1]-estimated[i-1].thrmv[0]*estimated[i-1].thrmv[1])/dl + (estimated[i+1].thrmv[1]-2*estimated[i].thrmv[1]+estimated[i-1].thrmv[1])/(dl*dl) );
		}
		dGdt[lpts] = exp(-2*estimated[lpts].l)*( PI*(-estimated[lpts-1].thrmv[0]*estimated[lpts-1].thrmv[1])/dl + (-2*estimated[lpts].thrmv[1]+estimated[lpts-1].thrmv[1])/(dl*dl) );

		/* Step G in time, for all length scales */
		for(i=1; i<=lpts; i++){
			estimated[i].thrmv[1] += dGdt[i];
		}

		/* Calculate state-variables at progressively larger length scales, up to coherence length, via Runge-Kutta method (derk) */
		for(i=1; i<=lpts; i++){
			// if(estimated[i].thrmv[0]>Kthresh) break;

			/* Calculate next estimated state (K-value) in time, for all length scales */
			estimated[i].thrmv[0] = derk1(IN FPRecRel, IN estimated[i-1], IN dl, IN N);
			estimated[i].tempv =  1.0-qTfrac;

			/* Print at most [datamax] (e.g., 1000) data points to screen and output file */
			if(tu<5 && (i%(lpts/datamax)==0||i<10)){   //  ((tu==1||tu==10||tu==100||tu==1000||tu==10000) && ((i+1)%(lpts/datamax)==0)){
				fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\n",
					1-estimated[i].tempv, estimated[i].tempv, estimated[i].l, estimated[i].thrmv[0], estimated[i].thrmv[1], estimated[i].thrmv[0]/K0);  // writes to a temporary file (not yet to outfile, without fflush)
				fflush(outfile);
				printf(         "%g\t%g\t%g\t%g\t%g\t%g\n",
					1-estimated[i].tempv, estimated[i].tempv, estimated[i].l, estimated[i].thrmv[0], estimated[i].thrmv[1], estimated[i].thrmv[0]/K0);
			}
		}
	}

	/* Close output file */
	fclose(outfile);
}



// vortex-pair (KT) theory in-equilibrium recursion relations
RETARRAY VPairRecRel(double l, double z[]){
	RETARRAY dzdl;

	dzdl.arr[0] = -B*z[1]*z[1]*z[0]*z[0];
	dzdl.arr[1] = (2-PI*z[0])*z[1];
	// dzdl.arr[2] = -2*PI*z[1]*exp(-2*l);

	return dzdl;
}

// vortex-pair (KT) theory out-of-equilibrium recursion relation
RETARRAY FPRecRel(double l, double z[]){
	RETARRAY dzdl;

	dzdl.arr[0] = -B*a04*z[0]*z[0]*exp(4*l)*z[1];
	dzdl.arr[1] = 0;
	// dzdl.arr[2] = -2*PI*z[1]*exp(-2*l);

	return dzdl;
}



/*  Program Notes  */
/*=========================================================================*/
/*

See further documentation (UCLAresearchNotes.pdf, UCLAresearchPrograms.pdf, derk1.c) for elaboration.


CONSTANTS
=========

PI is a familar mathematical constant.
PISQ is PI squared.
A is a dimension-dependent constant in the recursion relations (see dzdl.arr[0] in vltRecRel()).


PROGRAM PARAMETERS/INPUTS
=========================
lmax is...
dl is the step size for the increments in dimensionless length scale l...
K0c is the critical value of K0, which is the "bare K" or "K at smallest length scale a=a0".
   (K0c is K0 at T=Tc.)
Cc is the critical (lambda) vortex-loop-core energy (constant) (...for a loop of particular size?)
N: N+1 is the number of independent thermodynamic variables (?)...


VARIABLES
=========

The STATE structure was defined to group thermodynamic state information together:
state <-> {~temperature, ~length scale, superfluid ratio, fugacity, ~free energy}

Conceptually, there are:
* The real states; in reality, the system starts at an intial state and progresses through real states that can be measured approximately.
* The measured states; the data that are measured give a mathematical set of states that approximately describes the real states.
* The estimated states; using various numerical methods of differing accuracy, one can calculate states to estimate the measured states.
* The test states; in the 4th-order Runge-Kutta method, each step in calculating the next estimated state requires the calculation of three test states.
* The calculated states; the estimated and test states are calculated states (but the test states are not necessarily estimated states).

In this program (including derk1.c), there are "estimated", "calculated", and "test" states.

tempv = 1-T/Tc
    tempv is a "(dimensionless) temperature variable"
    T is the absolute temperature.
    Tc is the critical lambda temperature separating normal fluid from superfluid.

l is a dimensionless length scale related to a vortex-loop diameter a (l = ln(a/a0)).
    a is a vortex-loop diameter (which is larger than a0, the smallest such diameter).
    lmax is the largest value of l.
Given that D is the size of the system (the diameter of the spherical volume containing the helium), the program calculates to a maximum vortex diameter amax = D = a0*exp(lmax).

There are three thermodynamic variables (thermv[i] or z[i], i=1,2,3):
 thermv[0] and z[0] are K, the "dimensionless superfluid ratio" or "coupling constant".
 thermv[1] and z[1] are y, the fugacity.
 thermv[2] and z[2] are e, what I'll call the Helmholtz parameter.
In terms of other variables, these quantities are defined by the following relations:
  K = (hbar^2*a0*rhos)/(m^2*kB*T).
  y = exp{-mu/(kB*T)} = exp{-U/(kB*T)}.
  e = (A/kB*T)(a0^3/V), the Helmholtz free energy A divided by the product of Boltzmann constant kB with the absolute temperature T, multiplied by the volume a0^3 divided by the volume V of the system.  This is "the Helmholtz free energy per thermal energy per mini-system."  (The variable name "a" was already taken, so I chose "e".)

The RETARRAY structure was defined to allow the function vltRecRel to return an array of variables, which are derivatives of the thermodynamic variables.
 dzdl.arr[0] is dK/dl
 dzdl.arr[1] is dy/dl
 dzdl.arr[2] is de/dl
In derk1.c an array of RETARRAYs is defined
 DthrmvDl[i].arr[0] is dK/dl for the ith test STATE
 DthrmvDl[i].arr[1] is dy/dl for the ith test STATE
 DthrmvDl[i].arr[2] is de/dl for the ith test STATE


FUNCTIONS
=========

derk()
------
derk = "Differential Equation Runge Kutta method"
Dl is the increment in l.
n is the number of thermodynamic variables (n=3 given thermv[i], i=1,2,3).
  See derk1.c for more.

vltRecRel()
------------
vltRecRel = "vortex loop theory recursion relations"
Meaning of z[3] and zp[3]:
    z[0] = K  dzdl.arr[0] = dK/dl
    z[1] = y  dzdl.arr[1] = dy/dl
    z[2] = e  dzdl.arr[2] = de/dl


POTENTIAL FURTHER DEVELOPMENTS
==============================

Make pressure or volume a changable variable.

T P V
U H  A  G   fugacity y = exp(-u) = exp(-mu/kT); chemical potential mu
S S1 S2 S3
cV cP


*/























//			/* Bug Catching */ ///////////////////////////////////////////////////////////////////////////
//			if((tu==1||tu==10||tu==100||tu==1000||tu==10000) && ((i+1)%(lpts/datamax)==0 || i==1 || i==2)){
//				fprintf(outfile,"\t%i\t%g\t%g\n",
//					i,estimated[i].thrmv[1],dGdt[i] );
//				fflush(outfile);
//			}/////////////////////////////////////////////////////////////////////////////////////////////



//			/* Bug Catching */ ///////////////////////////////////////////////////////////////////////////
//			if((tu==1||tu==10||tu==100||tu==1000||tu==10000) && ((i+1)%(lpts/datamax)==0 || i==1 || i==2)){
//				fprintf(outfile,"\t\t%g\n",
//					estimated[i].thrmv[1] );
//				fflush(outfile);
//			}/////////////////////////////////////////////////////////////////////////////////////////////
