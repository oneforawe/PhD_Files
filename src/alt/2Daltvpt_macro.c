//  File Comments  //
//=========================================================================//

/* FILENAME: 2Daltvpt_macro.c */
/* VERSION: 1 (2011 Aug 03 - ...)
    In this program, I was going to try to use the equation derived in 4.2.2 of the Vortex Theories derivation paper (Method 3, the consistent method), but abandoned it since the equation has K0 in it (and Gary said that's bad).

    This is the 2D version of ring1.c and was formerly named 2DKT.c.  */
/* AUTHOR: Andrew Forrester <aforrester@ucla.edu> */
/* DESCRIPTION:
   * This C program calculates macroscopic (or system-scale) thermodynamic properties of a superliquid helium film at various temperatures, given a certain chosen system size lmax, and puts the results into one output file.
   * This program calculates macroscopic equilibrium thermodynamic (or quasi-thermo-static) properties of a 2D system of superfluid (a thin layer of liquid 4He), using the Kosterlitz-Thouless (KT) theory, which could be considered a "vortex pair theory" (vpt, as opposed to the 3D vortex loop theory, vlt), hence the name "2Daltvpt_macro.c".  By macroscopic, we mean that the quantities calculated are only recorded at the largest length scale, lmax (the size of the system), rather than recording all the renormalized values from the smallest scale (l=0 or a=a0) all the way up to lmax.
   * This program starts from a particular initial condition at the smallest length scale a0 and estimates the conditions at larger length scales for particular temperatures using differential recursion relations (obtained from the KT vortex-pair He-transition theory) and a 4th-order Runge-Kutta method.
   * The system is a thin film of liquid helium (4He) at atmospheric pressure, under or near the critical transition temperature Tc = Tkt separating normal fluid and superfluid phases.
   * The properties calculated are the quantities K, K/K0=sigma_s/sigma_0, y, and Gamma (not e?) as functions of temperature T (or Tfrac or tempv), which depend on the maximum length scale lmax (or amax) chosen.  See the Program Notes below the code for explanation of constants, variables, and functions.
*/
/* EXT FILES: none */
/* COMPILE NOTES:
   * To compile, type "g++ -lm 2Daltvpt_macro.c" without the quotes; then to run, type "./a.out".
   * NOTE: Be sure lsteps is divisible by ldatamax.
*/



//  Function Preparation  //
//=========================================================================//

// Standard routine header files //
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Constants //
const double PI   = 3.14159265358979323846;
const double PISQ = 9.86960440108935799230;
const double PICU = 31.00627668029981620634;
const double B    = 4.0*PICU;

// Program parameters/inputs //
const int    Op       = 3;         //  Option (1 -> hand-picked temp's; else temperature-spread formulas)
const double lmax     = 1.0;      //  1  //  maximum of the length scale l (potentially approximate)
const int    lsteps   = 1000;      //  1250, 5000, 100000 //  lsteps = lmax/dl = 10/0.0001 = 100,000  (make sure this is divisible by ldatamax)
const int    lpts     = lsteps+1;  //  from l=0 to l=lmax, inclusive (incl. the boundary-condition-enforcing pnt, there are lpts+1 "l" pnts)
const int    ldatamax = 100;       //  50,100 //  max number of data points printed (in bug check) per temperature examined
const double a0       = 1.0;       //  a0 in units of a0 (!)
const double a02      = a0*a0;     //  a0 to the second power
const double a04      = a02*a02;   //  a0 to the fourth power
const double K0c      = 0.747853;  //  critical (T=Tkt) bare superfluid ratio (or "coupling constant")



//  Function Prototypes  //
//=========================================================================//

// rk4((*f)(),*x,y[],h,n);
// EqRecRel(x,z,dzdx,n);
void rk4(void (*f)(double, double*, double*, unsigned int), double *x, double y[], double h, unsigned n, double Kzero);
void EqRecRel(double x, double z[2], double dzdx[2], unsigned n, double Knaught);



//  Function Definitions  //
//=========================================================================//

int main(){
	// Main function definitions //
	int    i,j, notdone=true, TempLabelOp;
	double dblsteps=lsteps, dl=lmax/dblsteps;
	double x, l[lpts];
	double Tfrac, tempv;
	double z[2], K[lpts], G[lpts];
	double K0, G0, y;
	FILE   *outfile;
	char   *filename;

	switch(Op){  //  TempLabelOp -> 0 or 1 for labeling the file with "T" or "Tv", respectively
		case 1:  TempLabelOp=0;  break;
		case 2:  TempLabelOp=1;  break;
		case 3:  TempLabelOp=0;  break;
		case 4:  TempLabelOp=0;  break;
		case 5:  TempLabelOp=1;  break;
	} // TempLabelOp, and the file label (T or Tv) is meant to indicate which quantity is best used for plotting the data T/Tc = T/Tkt or Tv = tempv = 1-T/Tc = 1-T/Tkt

	// Prepare output file, print identification and values //
	if(TempLabelOp==0)  asprintf(&filename, "2Daltvpt_macro_T_lmax%g_dl%g_Op%i.dat",  lmax,dl,Op);
	if(TempLabelOp==1)  asprintf(&filename, "2Daltvpt_macro_Tv_lmax%g_dl%g_Op%i.dat", lmax,dl,Op);
	outfile = fopen(filename,"w");  //  E.g., "2Daltvpt_constT.out"
	printf(     "\n\n# Filename: %s\n", filename);
	fprintf(outfile,"# Filename: %s\n", filename);
	fprintf(outfile,"# Source: 2Daltvpt_constT.c\n");
	fprintf(outfile,"# Source version: %s\n", "1 (2011 Aug 03 - ...)");
	fprintf(outfile,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g, B=%21.21g\n", PI,PISQ,PICU,B);
	fprintf(outfile,"# Parameter values: lmax=%g, lsteps=%i, (dl=%g), lpts=%i, ldatamax=%i, a0=%g, a02=%g, a04=%g, K0c=%g\n", lmax,lsteps,dl,lpts,ldatamax,a0,a02,a04,K0c);
	fprintf(outfile,"# Option values: Op=%i\n", Op);
	fprintf(outfile,"%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","y","G"); // what about e?
	printf(         "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","y","G");

	/* Perform calculations for various temperatures (calculations for each temp are independent) */
	i=0;
	while(notdone==true){
		switch(Op){
			case 1:	// A few hand-picked temperatures //
				switch(i){
					case 0:   Tfrac=1.00;  tempv=1.0-Tfrac;  i++;            break;
					case 1:   Tfrac=0.90;  tempv=1.0-Tfrac;  i++;            break;
					case 2:   Tfrac=0.80;  tempv=1.0-Tfrac;  i++;            break;
					case 3:   Tfrac=0.70;  tempv=1.0-Tfrac;  i++;            break;
					case 4:   Tfrac=0.60;  tempv=1.0-Tfrac;  i++;            break;
					case 5:   Tfrac=0.50;  tempv=1.0-Tfrac;  i++;            break;
					case 6:   Tfrac=0.40;  tempv=1.0-Tfrac;  i++;            break;
					case 7:   Tfrac=0.30;  tempv=1.0-Tfrac;  i++;            break;
					case 8:   Tfrac=0.20;  tempv=1.0-Tfrac;  i++;            break;
					case 9:   Tfrac=0.10;  tempv=1.0-Tfrac;  notdone=false;  break;
				}
				break;
			case 2:	// Decrement below Tc (staying very close to Tc) //
				// from T/Tc = 1-1e-8 to T/Tc = 0.992
				tempv = 0.0+pow(10,-8.0+0.1*i);  Tfrac = 1.0-tempv;
				if(i==59) notdone=false;
				i++;  break;
			case 3:	// Increment from near zero to above Tc //
				// from T/Tc = 0.206 to T/Tc = 0.9990;  then from T/Tc = 1.0012 to T/Tc = 2
				if(i<30){ tempv = 0.0+pow(10,-0.1-0.1*i);      Tfrac = 1.0-tempv; }
				else    { tempv = 0.0-pow(10,-5.0+0.1*(i-9));  Tfrac = 1.0-tempv; }
				if(i==59) notdone=false;
				i++;  break;
			case 4:	// Increment from near zero to just below Tc //
				// from T/Tc = 0.206 to T/Tc = 0.9991
				tempv = 0.0+pow(10,-0.1-0.05*i);  Tfrac = 1.0-tempv;
				if(i==59) notdone=false;
				i++;  break;
			case 5:	// Decrement below Tc (staying extremely close to Tc) //
				// T/Tc = 1-1e-20  to  T/Tc ~ 1-1e-16
				tempv = 0.0+pow(10,-20.0+0.1*i);  Tfrac = 1.0-tempv;
				// earlier: T/Tc = 1-1e-12  to  T/Tc ~ 1-1e-6
				if(i==59) notdone=false;
				i++;  break;
		}

		// Initialize quantities at smallest length-scale (l=0) //
		x = l[0] = 0.0;
		z[0] = K[0] = K0 = K0c/(1.0-tempv);
		z[1] = G[0] = G0 = exp(-4.0*l[0]-PISQ*K0)/a04;

		// Calculate state-variables at progressively larger length scales, up to coherence length, via Runge-Kutta method (rk4) //
		for(j=1; j<lpts; j++){
			// Step out in length scale l, calculating next l, K, and G //
			rk4(EqRecRel, &x, z, dl, 2, K0);  // equil: rk4 2D (K and G), EqRecRel
			l[j] = l[j-1] + dl;
			K[j] = z[0];
			G[j] = z[1];
			y    = a02*exp(2.0*l[j])*sqrt(G[j]);
			// Bug Check: Print at most [ldatamax] (e.g., 50) data points to screen and output file //
			//if(j%(lsteps/ldatamax)==0){
			//	printf(         "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", Tfrac, tempv, l[j], K[j], K[j]/K0, y, G[j]);
			//	fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\t%g\n", Tfrac, tempv, l[j], K[j], K[j]/K0, y, G[j]);
			//}
		}

		/* Record data values */
		printf(         "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", Tfrac, tempv, l[lpts-1], K[lpts-1], K[lpts-1]/K0, y, G[lpts-1]);
		fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\t%g\n", Tfrac, tempv, l[lpts-1], K[lpts-1], K[lpts-1]/K0, y, G[lpts-1]);
		fflush(outfile);
	}

	fclose(outfile);

}



// Runge-Kutta 4th order method //
void rk4(void (*f)(double, double*, double*, unsigned), double *x, double y[], double h, unsigned n, double Kzero){
	int i,j;
	double *k[4],*s,temp;

	for(i=0;i<4;i++)  k[i] = (double *)calloc(n,sizeof(double));
	s = (double *)calloc(n,sizeof(double));
	f(*x,y,k[0],n,Kzero);
	for(i=1;i<4;i++){
		temp = h*((i+1)/2)/2;
		// using integer division trick (times double): ((i+1)/2)/2 = 0.50, 0.50, 1.00 instead of 0.75, 1.00, 1.25
		for(j=0;j<n;j++)  s[j] = y[j]+k[i-1][j]*temp;
		f(*x+h*(i/2)/2,s,k[i],n,Kzero);
		// using integer division trick (times double): (i/2)/2 = 0.00, 0.50, 0.50 instead of 0.50, 0.75, 1.00
	}
	for(j=0;j<n;j++)  y[j] += h*(k[0][j]+2.0*(k[1][j]+k[2][j])+k[3][j])/6.0;
	*x += h;
	for(i=0;i<4;i++)  free((char *) k[i]);
	free((char *)s);
}



// equilibrium (K,y) recursion relations //
void EqRecRel(double x, double z[2], double dzdx[2], unsigned n, double Knaught){
	dzdx[0] = -B*z[1]*z[1]*z[0]*z[0]*z[0]*z[0]/Knaught/Knaught/(1.0+B*z[1]*z[1]*z[0]*z[0]*z[0]/Knaught/Knaught);  // Wait! this eqn has K0 in it... that's probably a bad thing...
	dzdx[1] = z[1]*(2.0-PI*z[0]);
}



//  Program Notes  //
//=========================================================================//
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

There are three thermodynamic variables (thermv[i], i=1,2,3):
thermv[0] is K, the "coupling constant", or "dimensionless superfluidity ratio".
thermv[1] is y, the fugacity (e^{-mu/kT}).
thermv[2] is F, the (Helmholtz) free energy of the system.


FUNCTIONS
=========

derk()
------
Dl is the increment in l.
n is the number of thermodynamic variables (n=3 given thermv[i], i=1,2,3).
  See derk1.c for more.

ringrecrel()
------------
Meaning of z[3] and zp[3]:
    z[0] = K    dzdl.arr[0] = dK/dl
    z[1] = y    dzdl.arr[1] = dy/dl
    z[2] = F    dzdl.arr[2] = dF/dl

*/
