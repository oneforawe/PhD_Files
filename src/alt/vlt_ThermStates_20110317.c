//  File Comments  //
//=========================================================================//

/* FILENAME: vlt_ThermStates.c */
/* VERSION: 5 (2011 Mar 01 - 2011 Mar 17)
/* AUTHOR: Andrew Forrester <aforrester@ucla.edu> */
/* DESCRIPTION:
   * This C program calculates (length-scale-dependent) thermodynamic states of a system of liquid helium-4 (4He) under conditions specified by user inputs.
   * The system is a spherical region (of variable diameter greater than D) of liquid 4He at some constant pressure, near or below the critical transition temperature Tc (which depends on pressure) separating normal fluid and superfluid (aka the lambda temperature).
   * A state is merely a set of values of certain quantities that correspond to temperature T, length scale a, dimensionless superfluid ratio K, dimensionless fugacity y, and dimensionless Helmholtz energy parameter e.  The "bare" superfluid ratio K0 at the smallest relevant length scale a0 (the diameter of the smallest vortex loops) and the renormalized superfluid ratio Kr (and Kr/K0) are also of interest.  See the "Program Notes" below the code for explanation of these quantities and the constants, variables, and functions in this program.
   * The actual values of the set of temperatures {Tn} and the length scale a are irrelevant; the calculations are made using unitless variables (tempv and l) scaled by Tc (which depends on pressure) and a0 (which depends on temperature and pressure).  The actual value of the pressure is also irrelevant; the user may choose a value of Cc without knowing what pressure it corresponds to.
   * The program uses Gary A. Williams's vortex-loop theory of the superfluid phase transition (i.e., its differential recursion relations) and a 4th-order Runge-Kutta integration method to calculate (length-scale-dependent) thermodynamic quantities (essentially K0,K,Kr,y,e) of the system at a set of temperatures {Tn}, at various length scales (l, no greater than a length scale corresponding to the choice of lmax), and at some pressure (corresponding to the choice of Cc, the critical dimensionless vortex-loop-core energy).
   * For a set of temperatures {Tn}, the program starts from some initial state (l=0,K0=K0(Cc,Tn),K=K0,Kr=K0,y=y(K0,Cc),e=0) (depending on each temperature Tn and the choice of Cc) at a smallest length scale (l=0 and a=a0), it integrates to larger length scales (incrementing by Dl using the Runge-Kutta technique) up to some limiting length scale (see while loop, where l<lmax and K<5), and it estimates and records the state at that final length scale (l=lf,K0=K0(Cc,Tn),K=Kf,Kr=Kr(K,l),y=yf,e=ef).
   * The various temperatures are determined by "Options" (AA, A, B, and C) in the program, where one can choose (by manually editing this file) sets of temperatures below, above, or very near Tc.  One could also add more options to choose different sets of temperatures.

   * In the basename "vlt_ThermStates", "vlt" refers to Williams's Vortex Loop Theory (of the superfluid phase transition) and "ThermStates" refers to the THERModynamic STATES described by the calculated thermodynamic quantities.  (This file was previously called ring3.c.)

   Inputs:
   * N       - You can choose N=2 (if you don't care to calculate e) or N=3.
   * Kthresh -
   * dl      -
   * lmax    - Choosing 100 corresponds to a system with diameter no greater than D = a0*exp(lmax) ~ 6.7e33 meters (given that a0 ~ 2.5 angstroms at low pressures), allowing for a huge system and calculations at enormous length scales.
   * From vlt_K0cFind.c, you can get a set of {Cc, K0c} data to put in the 
   * Select (uncomment) "Option" (AA, A, B, or C) to choose a set of temperatures.

   Output:
   * This program returns text files of data in a two dimensional arrays.

   Questions:
   1. Should we say that the system is of variable volume or should we say that the system is of infinite volume?  (D would then just be a limiting diameter for the calculations.)  If we are assuming that the system is of fixed, finite volume (diameter D), then when you select the temperature of the system, doesn't that also determine the pressure (and thus, wouldn't that determine a certain Cc for each temperature)?
   2. Why do you use
      estimated.thrmv[0]<5.0
   in this program but
      estimated.thrmv[0]<6.0
   in vlt_HeatCap_P.c?
   3. Why don't you include a finite-size option in vlt_HeatCap_P.c
      while(estimated.thrmv[0]<6)
   like you do in this file:
      while(estimated.l<lmax && estimated.thrmv[0]<5.0) ?
*/
/* EXT FILES: none */
/* COMPILE NOTES:
   * To compile, type "g++ -lm vlt_ThermStates.c" without the quotes; then to run, type "./a.out".
*/



//  Function Preparation  //
//=========================================================================//

// Standard routine header files //
#include <stdio.h>
#include <math.h>

// Constants and labels //
#define IN           // input label
#define OUT          // output label
#define INOUT        // input/output label
const double PI    = 3.14159265358979323846;
const double PISQ  = 9.86960440108935861883;
const double A     = 4.0*PI*PISQ/3.0;
const double THETA = 0.6;

// Program parameters/inputs //
const int    N       = 3;
const double Kthresh = 5.0;
const double dl      = 0.0001;
const int    lmax    = 100;    // 100 corresponds to largest length scale D = a0*exp(lmax) ~ 6.8e33 meters (huge diameter system)

// Data type definitions //
typedef struct {
	double l;
	double tempv;
	double thrmv[N];
} STATE;

typedef struct {
	double arr[N];
} RETARRAY;



//  Function Prototypes  //
//=========================================================================//

// vltRecRel(l,z)
// rk4(func(),calculated,Dl,n)
RETARRAY vltRecRel(double l, double z[], double C);
STATE rk4(RETARRAY (*func)(double l, double z[], double C), STATE calculated, double Dl, unsigned n, double c);



//  Function Definitions  //
//=========================================================================//

main(){
	// Main function definitions //
	int s,i,j=0;
	double Cc,K0c,K0,Kr;
	STATE estimated;
	FILE *outfile;
	char *filename;

	// Progress through 22 sets of values of Cc and K0c, where a data file will be produced for each set //
	for(s=1; s<=25; s++){
		switch(s){
			// from vlt_K0cFind.c using vltRecRel with Cc: //
			case 1:		Cc = 0.29;	K0c = 0.590412770819029;	break;
			case 2:		Cc = 0.285;	K0c = 0.598791680907394;	break;
			case 3:		Cc = 0.28;	K0c = 0.60744048117457;		break;
			case 4:		Cc = 0.275;	K0c = 0.616373127790686;	break;
			case 5:		Cc = 0.27;	K0c = 0.625604581604229;	break;
			case 6:		Cc = 0.265;	K0c = 0.635150901567938;	break;
			case 7:		Cc = 0.26;	K0c = 0.645029348886667;	break;
			case 8:		Cc = 0.255;	K0c = 0.655258503360094;	break;
			case 9:		Cc = 0.25;	K0c = 0.665858393630086;	break;
			case 10:	Cc = 0.245;	K0c = 0.676850643323299;	break;
			case 11:	Cc = 0.24;	K0c = 0.688258635413654;	break;
			case 12:	Cc = 0.235;	K0c = 0.700107697527901;	break;
			case 13:	Cc = 0.23;	K0c = 0.712425311394914;	break;
			case 14:	Cc = 0.225;	K0c = 0.72524135021341;		break;
			case 15:	Cc = 0.22;	K0c = 0.738588348405702;	break;
			case 16:	Cc = 0.215;	K0c = 0.752501809064936;	break;
			case 17:	Cc = 0.21;	K0c = 0.767020555425591;	break;
			case 18:	Cc = 0.205;	K0c = 0.782187133937134;	break;
			case 19:	Cc = 0.2;	K0c = 0.79804827805642;		break;
			case 20:	Cc = 0.195;	K0c = 0.814655443770591;	break;
			case 21:	Cc = 0.19;	K0c = 0.83206543021527;		break;
			case 22:	Cc = 0.185;	K0c = 0.85034110168895;		break;
			case 23:	Cc = 0.18;	K0c = 0.869552231049103;	break;
			case 24:	Cc = 0.175;	K0c = 0.889776489127531;	break;
			case 25:	Cc = 0.17;	K0c = 0.911100610712701;	break;
			/* from vlt_K0cFind.c using vltRecRel with Cc: //
			case 1:		Cc = 1.20;	K0c = 0.181368778098505;	break;
			case 2:		Cc = 1.10;	K0c = 0.195539960611704;	break;
			case 3:		Cc = 1.06;	K0c = 0.201872397682003;	break;
			case 4:		Cc = 1.05;	K0c = 0.203522562154756;	break;
			case 5:		Cc = 1.04;	K0c = 0.205200971992232;	break;
			case 6:		Cc = 1.03;	K0c = 0.206908384317726;	break;
			case 7:		Cc = 1.02;	K0c = 0.208645584189656;	break;
			case 8:		Cc = 1.01;	K0c = 0.210413385921907;	break;
			case 9:		Cc = 1.00;	K0c = 0.212212634480595;	break;
			case 10:	Cc = 0.99;	K0c = 0.214044206962472;	break;
			case 11:	Cc = 0.98;	K0c = 0.215909014160657;	break;
			case 12:	Cc = 0.97;	K0c = 0.217808002223811;	break;
			case 13:	Cc = 0.90;	K0c = 0.232145749019357;	break;
			case 14:	Cc = 0.80;	K0c = 0.25645126154457;		break;
			case 15:	Cc = 0.70;	K0c = 0.28683227513254;		break;
			case 16:	Cc = 0.60;	K0c = 0.326048538761259;	break;
			case 17:	Cc = 0.55;	K0c = 0.350345580223672;	break;
			case 18:	Cc = 0.50;	K0c = 0.378915731974668;	break;
			case 19:	Cc = 0.40;	K0c = 0.454738952186452;	break;
			case 20:	Cc = 0.30;	K0c = 0.574413324380367;	break;
			case 21:	Cc = 0.20;	K0c = 0.79804827805642;		break;
			case 22:	Cc = 0.10;	K0c = 1.41499685107337;		break;
			*/
			/* from vlt_K0cFind.c using vltRecRel without Cc:
			case 1:		Cc = 1.20;	K0c = 0.282467106542233;	break;
			case 2:		Cc = 1.10;	K0c = 0.298805159967420;	break;
			case 3:		Cc = 1.06;	K0c = 0.306075989642654;	break;
			case 4:		Cc = 1.05;	K0c = 0.307968028765089;	break;
			case 5:		Cc = 1.04;	K0c = 0.309891387590691;	break;
			case 6:		Cc = 1.03;	K0c = 0.311846911030558;	break;
			case 7:		Cc = 1.02;	K0c = 0.313835475596892;	break;
			case 8:		Cc = 1.01;	K0c = 0.315857990917745;	break;
			case 9:		Cc = 1.00;	K0c = 0.317915401340550;	break;
			case 10:	Cc = 0.99;	K0c = 0.320008687630639;	break;
			case 11:	Cc = 0.98;	K0c = 0.322138868771385;	break;
			case 12:	Cc = 0.97;	K0c = 0.324307003873199;	break;
			case 13:	Cc = 0.90;	K0c = 0.340645049254360;	break;
			case 14:	Cc = 0.80;	K0c = 0.368242208355569;	break;
			case 15:	Cc = 0.70;	K0c = 0.402631447770970;	break;
			case 16:	Cc = 0.60;	K0c = 0.446957413953460;	break;
			case 17:	Cc = 0.55;	K0c = 0.474427344215761;	break;
			case 18:	Cc = 0.50;	K0c = 0.506764237544466;	break;
			case 19:	Cc = 0.40;	K0c = 0.592903866733522;	break;
			case 20:	Cc = 0.30;	K0c = 0.730151410553482;	break;
			case 21:	Cc = 0.20;	K0c = 0.991358025249279;	break;
			case 22:	Cc = 0.10;	K0c = 1.738147692720230;	break;
			*/
		}

		// Prepare output file, print identification, values, and headings for data (to screen, too) //
		asprintf(&filename, "vlt_ThermStates_lmax%i_A_Cc%4.3f.dat", lmax, Cc);
		outfile = fopen(filename,"w");  //  E.g., "vlt_ThermStates.out" or "vlt_ThermStates_lmax100_A_Cc0.00.dat"
		fprintf(outfile,"# Filename: %s\n", filename);
		fprintf(outfile,"# Source: vlt_ThermStates.c\n");
		fprintf(outfile,"# Source version: %s\n", "5 (2011 Mar 01 - ...)");
		fprintf(outfile,"# Constant values: PI=%21.21g, PISQ=%21.21g, A=%21.21g\n", PI,PISQ,A);
		fprintf(outfile,"# Parameter values: N=%i, Kthresh=%g, dl=%g, lmax=%i\n", N,Kthresh,dl,lmax);
		fprintf(outfile,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","y","e","Kr","Kr/K0");
		printf(         "\n%s%g\t%s%g%s\n\n", " ---  Cc = ",Cc,"K0c = ",K0c,"  ---");
		printf(         "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","y","e","Kr","Kr/K0");

		// Loop through different temperatures (calculations for each temperature are independent) //
		for(i=0; i<60; i++){
			// (Option AA) Set temperature and comment-out for-loop with index i //
			//estimated.tempv = 0.0;
			// (Option A) Decrement temperature (further below Tc) //
			estimated.tempv = 0.0+pow(10,-8.0+0.1*i);
			// (Option B) Increment temperature (near zero to above Tc) //
			//if(i<30) estimated.tempv = 0.0+pow(10,-0.1-0.1*i);
			//else     estimated.tempv = 0.0-pow(10,-5.0+0.1*(i-9)); */
			// (Option C) Increment temperature (near zero to just below Tc) //
			//estimated.tempv = 0.0+pow(10,-0.1-0.05*i);

			// Initialize length scale and energy values //
			estimated.l = 0.0;
			estimated.thrmv[2] = 0.0;

			// Initialize temperature-dependent values at initial length scale //
			K0 = K0c/(1.0-estimated.tempv);
			estimated.thrmv[0] = K0;
			estimated.thrmv[1] = exp(-PISQ*K0*Cc);  //  1.0/(exp(PISQ*K0*Cc)-1.0);

			// Calculate state-variables at progressively larger length scales, up to coherence length, via Runge-Kutta method (rk4) //
			while(estimated.l<lmax && estimated.thrmv[0]<Kthresh){
				// Calculate next estimated state, incrementing length scale //
				estimated = rk4(IN vltRecRel, INOUT estimated, IN dl, IN N, IN Cc);
				/* Print intermediate numbers to screen (for bug checking)
				j++;
				if(j%2000==0){
					Kr = estimated.thrmv[0]*exp(-estimated.l);
					//fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
					//	1-estimated.tempv, estimated.tempv, estimated.l, estimated.thrmv[0], estimated.thrmv[1], estimated.thrmv[2], Kr, Kr/K0);  // writes to a temporary file (not yet to outfile, see fflush below)
					//fflush(outfile);  // forces write-to-file (outfile) from the temporary file (see fprintf above), so if an error interrupts the execution of the code at least some data will be recorded (for bug checking)
					printf(         "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
						1-estimated.tempv, estimated.tempv, estimated.l, estimated.thrmv[0], estimated.thrmv[1], estimated.thrmv[2], Kr, Kr/K0);
				}
				if(j>100000) j=0; */
			}

			// Update values that depend on state-variables, if not using "Print out intermediate numbers" block of code //
			Kr = estimated.thrmv[0]*exp(-estimated.l);

			// Print data to screen and output file //
			fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
				1-estimated.tempv, estimated.tempv, estimated.l, estimated.thrmv[0], estimated.thrmv[1], estimated.thrmv[2], Kr, Kr/K0);  // writes to a temporary file (not yet to outfile, see fflush above)
			printf(         "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
				1-estimated.tempv, estimated.tempv, estimated.l, estimated.thrmv[0], estimated.thrmv[1], estimated.thrmv[2], Kr, Kr/K0);
		}

		// Close output file //
		fclose(outfile);
	}
}



// vortex loop theory recursion relations
RETARRAY vltRecRel(double l, double z[], double C){
	RETARRAY dzdl;

	dzdl.arr[0] = z[0]-A*z[0]*z[0]*z[1];
	dzdl.arr[1] = z[1]*( 6.0 - PISQ*z[0]*(1.0+C-THETA*log(z[0])) );
	// dzdl.arr[2] = -PI*z[1]*exp(-3.0*l);

	return dzdl;
}



// Runge-Kutta 4th order method //
STATE rk4(RETARRAY (*func)(double l, double z[], double C), STATE calculated, double Dl, unsigned n, double c){
	int i,j;
	double DL;
	RETARRAY DthrmvDl[4];
	STATE test[3];

	// Calculate the derivatives at the input calculated state and 3n test states //
	DthrmvDl[0] = func(IN calculated.l, IN calculated.thrmv, IN c);
	for(i=0; i<3; i++){
		DL = Dl*((i+2)/2)/2;                        // See "Note" below.
		for(j=0; j<n; j++)
			test[i].thrmv[j] = calculated.thrmv[j] + DthrmvDl[i].arr[j]*DL;
		test[i].l = calculated.l + Dl*((i+1)/2)/2;  // See "Note" below.
		DthrmvDl[i+1] = func(IN test[i].l, IN test[i].thrmv, IN c);
	}
	// Note: This algorithm, using '((...)/2)/2', creates a desired rounding "error" that produces the sequence 0.0, 0.5, 0.5, 1.0, instead of 0.25, 0.50, 0.75, 1.00.

	// Update the calculated state, using the derivatives calculated above //
	for(j=0; j<n; j++) calculated.thrmv[j] += Dl*(DthrmvDl[0].arr[j]+2*(DthrmvDl[1].arr[j]+DthrmvDl[2].arr[j])+DthrmvDl[3].arr[j])/6;
	calculated.l += Dl;

	// Return the new calculated state //
	return calculated;
}



//  Program Notes  //
//=========================================================================//
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

rk4()
------
rk4 = "Runge Kutta method, 4th order"
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
