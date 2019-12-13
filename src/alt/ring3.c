/* filename: ring3.c */
/* description: This program calculates thermodynamic properties of a system of liquid helium-4 at different temperatures and length scales.  It starts from a particular initial condition at the smallest length scale a0 and estimates the conditions at larger length scales for particular temperatures using differential recursion relations (obtained from a vortex-loop He-transition theory) and a 4th-order Runge-Kutta method.
       The system is a spherical volume of diameter D of liquid helium (4He) at atmospheric pressure, slightly under the critical transition temperature Tc (2.1768 K) separating normal fluid and superfluid.
       The properties calculated are the quantities K, y, and F (and Kr/K0 and a0/Kr) as functions of temperature T (or tempv) and length scale l (or a).  See below the code for explanation of constants, variables, and functions.
       The filename ring refers to the vortex rings (or loops) in this theory.
       (The program ring1.c is my first working nontrivial variation on the ring.c program that Dr. Gary Williams wrote.)
 */
/* ext files: derk1.c */
/* author: Andrew Forrester <aforrester@ucla.edu> */

// NOTE: Be sure that N=3 in derk1.c when compiling.


#include <stdio.h>
#include <math.h>

#define IN
#define OUT
#define INOUT
#define PI    3.14159265358979323846
#define PISQ  9.86960440108935861883
#define A     4.0*PI*PI*PI/3.0
#define kB    1.380650424e-23
#define hbar  1.05457162853e-34
#define m     6.646476186672e-27
#define a0    2.53
#define a0v   2.53e-10  // a0 value in meters
#define Tc    2.1768
#define K0c   0.3091469984216758
#define Cc    1.03      // 0.73109405147652 ?
#define N     3

/*
Cc	K0c (from K0cFind1.c)
1.20	0.2805443089651988
1.10	0.2964689389987769
1.03	0.3091469984216758
1.00	0.3150360068926090
0.90	0.3370353710663324
0.80	0.3636228655823345
0.70	0.3965669728967951
0.60	0.4387310072100319
0.55	0.4646970470439841
0.50	0.4951082899772983
0.40	0.5753434296585114
0.30	0.7010909420810494
0.20	0.9346549781449027
0.10	1.576187055119766
0.00	5.105045166053440
*/
/*
// Brooks and Donnelly fit:
#define D1    17.41647
#define D2    -60.48823
#define D3    -0.5307478
#define D4    1.817261e4
#define D5    -1.351398e7
#define D6    1.621499e8
#define D7    -5.062661e8
*/
// My fit for :
#define D1    10.1169430121868
#define D2    -9.88264499497637
#define D3    -0.664377974941481
#define D4    18200.0901963049
#define D5    -13499999.9767995
#define D6    162000000.003397
#define D7    -505999999.999502

#define r0    0.14514
#define r1    -0.000103362
#define r2    0.000340383
#define r3    -0.000450864
#define r4    0.000182967
#define rs0   145.12
#define rs1   0.0670526
#define rs2   -0.0556804
#define rs3   -0.111639
#define rs4   0.0687962

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
//   ringrecrel(l,z)

extern STATE derk(RETARRAY (*func)(double l, double z[]), STATE calculated, double Dl, unsigned n);

RETARRAY ringrecrel(double l, double z[]);



/*  Function definitions  */
/*=========================================================================*/

main(){
	int i,j=0;
	int Dexp=100;  // 100 corresponds to D = 6.8e33 meters (huge system)
	double T,rho,rhos,x,Delok,C,K0,Kr,a;
	STATE estimated;
	FILE *outfile;

	/* Print headings for data to screen and output file */
	outfile = fopen("ring3_Dexp100_A_beta_Cc1.20.dat","w");  //  "ring3.out"
	fprintf(outfile,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
		"T/Tc","1-T/Tc","l","a","K","y","F","Kr","Kr/K0","a0/Kr","Del/kB","C");
	printf(         "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
		"T/Tc","1-T/Tc","l","a","K","y","F","Kr","Kr/K0","a0/Kr","Del/kB","C");

	/* Loop through different temperatures (calculations for each temp are independent) */
	for(i=0; i<60; i++){
		/* (Option A) Decrement temperature (further below Tc) */
		estimated.tempv = 0.0+pow(10,-8.0+0.1*i);
		/* (Option B) Increment temperature (near zero to above Tc)
		if(i<30) estimated.tempv = 0.0+pow(10,-0.1-0.1*i);
		else     estimated.tempv = 0.0-pow(10,-5.0+0.1*(i-9)); */
		/* (Option C) Increment temperature (near zero to just below Tc)
		estimated.tempv = 0.0+pow(10,-0.1-0.05*i); */

		/* Initialize lengthscale and energy values */
		estimated.l = 0.0;
		estimated.thrmv[2] = 0.0;

		/* Initialize temperature-dependent values at initial lengthscale */
		T = (1-estimated.tempv)*Tc;
		rho = r0 + r1*T + r2*T*T + r3*T*T*T + r4*T*T*T*T;
		rhos = rs0 + rs1*T + rs2*T*T + rs3*T*T*T + rs4*T*T*T*T;
		x = -(D1+D2*rho)/T;
		  // Option alpha: simple factor multiplication
		//Delok = (6.8/5.33150354476608)*(D1 + D2*rho + D3*exp(x)*T/rho + D4*exp(2*x)*rho + (D5+D6*rho+D7*rho*rho)*exp(3*x));
		  // Option beta: complicated factor multiplication
		Delok = ((1.31086e-15)*exp(34.622*(1-estimated.tempv)-1.64327)+0.99999)*(D1 + D2*rho + D3*exp(x)*T/rho + D4*exp(2*x)*rho + (D5+D6*rho+D7*rho*rho)*exp(3*x));
		C = m*m*kB*Delok/(PISQ*hbar*hbar*a0v*rhos);  // (we set =1.3, value at 1.5K) (new comment: what? why 1.5K instead of Tc=2.1768K?)
/*		rho = r0+r1*T+r2*pow(T,2)+r3*pow(T,3)+r4*pow(T,4);
		rhos = rs0+rs1*T+rs2*pow(T,2)+rs3*pow(T,3)+rs4*pow(T,4);
		x = -(D1+D2*rho)/T;
		Delok = D1 + D2*rho + D3*exp(x)*T/rho + D4*exp(2*x)*rho + (D5+D6*rho+D7*pow(rho,2))*exp(3*x);
		C = pow(m,2)*kB*Delok/(PISQ*pow(hbar,2)*a0v*rhos); */
		K0 = K0c/(1.0-estimated.tempv);
		estimated.thrmv[0] = K0;
		estimated.thrmv[1] = exp(-PISQ*K0*Cc);

		/* Progress to larger lengthscales via Runge-Kutta method */
		while(estimated.l<Dexp && estimated.thrmv[0]<5.0){
			/* Calculate next estimated state, incrementing lengthscale */
			estimated = derk(IN ringrecrel, INOUT estimated, IN 0.0001, IN N);  // Dl=0.0001, n=2 or 3
			/* Print out intermediate numbers
			a = a0*exp(estimated.l);
			Kr = estimated.thrmv[0]*a0/a;
			j++;
			if(j%2000==0){
				fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
					1-estimated.tempv, estimated.tempv, estimated.l, a, estimated.thrmv[0], estimated.thrmv[1], estimated.thrmv[2], Kr, Kr/K0, a0/Kr);  // writes to a temporary file (not yet to outfile, see fflush below)
				printf(         "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
					1-estimated.tempv, estimated.tempv, estimated.l, a, estimated.thrmv[0], estimated.thrmv[1], estimated.thrmv[2], Kr, Kr/K0, a0/Kr);
			}
			if(j>100000)  j=0; */
		}

		/* Update dependent values, if not using "Print out intermediate numbers" block of code */
		a = a0*exp(estimated.l);
		Kr = estimated.thrmv[0]*a0/a;

		/* Record data values */
		fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
			1-estimated.tempv, estimated.tempv, estimated.l, a, estimated.thrmv[0], estimated.thrmv[1], estimated.thrmv[2], Kr, Kr/K0, a0/Kr,Delok,C);  // writes to a temporary file (not yet to outfile, see fflush below)
		fflush(outfile);  // forces write-to-file (outfile) from the temporary file (see fprintf above), so if an error interrupts the execution of the code at least some data will be recorded (for bug tracking)
		printf(         "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
			1-estimated.tempv, estimated.tempv, estimated.l, a, estimated.thrmv[0], estimated.thrmv[1], estimated.thrmv[2], Kr, Kr/K0, a0/Kr,Delok,C);
	}

	fclose(outfile);
}



// vortex "ring" recursion relations
RETARRAY ringrecrel(double l, double z[]){
	RETARRAY dzdl;

	dzdl.arr[0] = z[0]-A*z[1]*z[0]*z[0];
	dzdl.arr[1] = z[1]*(6.0-PISQ*z[0]*(1.0-0.6*log(z[0])));
	dzdl.arr[2] = -PI*z[1]*exp(-3.0*l);

	return dzdl;
}



/*

See further documentation (UCLAresearchNotes.pdf, UCLAresearchPrograms.pdf, derk1.c) for elaboration.


CONSTANTS
=========

PI is a familar mathematical constant.
PISQ is PI squared.
a0 is smallest vortex-loop diameter, in angstroms (10^{-10} meters).
K0c is the critical value of K0, which is the "bare K" or "K at smallest length scale a=a0".
   (K0c is K0 at T=Tc.)
C is the vortex-loop-core energy (constant) (...for a loop of particular size?)
A is a dimension-dependent constant in the recursion relations (see dzdl.arr[0] in ringrecrel()).


VARIABLES
=========

The STATE structure was defined to group thermodynamic state information together:
state <-> {temperature, length scale, superfluidity ratio, fugacity, free energy}

Conceptually, there are:
* The real states; in reality, the system starts at an intial state and progresses through real states that can be measured approximately.
* The measured states; the data that are measured give a mathematical set of states that approximately describes the real states.
* The estimated states; using various numerical methods of differing accuracy, one can calculate states to estimate the measured states.
* The test states; in the 4th-order Runge-Kutta method, each step in calculating the next estimated state requires the calculation of three test states.
* The calculated states; the estimated and test states are calculated states.

In this program (including derk1.c), there are "estimated" "calculated" and "test" states.

tempv = 1-T/Tc
    tempv is a "(dimensionless) temperature variable"
    T is the absolute temperature.
    Tc is the critical temperature separating normal fluid from superfluid.

l is a dimensionless length scale related to a vortex ring diameter (l = ln(a/a0)).
    a is a vortex-loop diameter (which is larger than a0, the smallest such diameter).
    Dexp is the largest value of l.
Given that D is the size of the system (the diameter of the spherical volume containing the helium), the program calculates to a maximum vortex diameter a_max = D.  Thus, we have l_max = ln(a_max/a0) = ln(D/a0) := Dexp, and D = a0*e^Dexp.  (That last equation shows why the name "Dexp" makes sense.)

There are three thermodynamic variables (thermv[i] or z[i], i=1,2,3):
 thermv[0] and z[0] are K, the "coupling constant", or "dimensionless superfluidity ratio".
 thermv[1] and z[1] are y, the fugacity (e^{-mu/kT}).
 thermv[2] and z[2] are F, the (Helmholtz) free energy of the system.

The RETARRAY structure was defined to allow the function ringrecrel to return an array of variables, which are derivatives of the thermodynamic variables.
 dzdl.arr[0] is dK/dl
 dzdl.arr[1] is dy/dl
 dzdl.arr[2] is dF/dl
In derk1.c an array of RETARRAYs is defined
 DthrmvDl[i].arr[0] is dK/dl for the ith test STATE
 DthrmvDl[i].arr[1] is dy/dl for the ith test STATE
 DthrmvDl[i].arr[2] is dF/dl for the ith test STATE


FUNCTIONS
=========

derk()
------
derk = "Differential Equation Runge Kutta method"
Dl is the increment in l.
n is the number of thermodynamic variables (n=3 given thermv[i], i=1,2,3).
  See derk1.c for more.

ringrecrel()
------------
ringrecrel = "vortex 'ring' recursion relations"
Meaning of z[3] and zp[3]:
    z[0] = K    dzdl.arr[0] = dK/dl
    z[1] = y    dzdl.arr[1] = dy/dl
    z[2] = F    dzdl.arr[2] = dF/dl

*/
