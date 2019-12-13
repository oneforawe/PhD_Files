//  File Comments  //
//=========================================================================//

/* FILENAME: 2Dvpt_K0cFind.c */
/* VERSION: 1 (2011 Dec 14 - 2012 Jul 31)
            First version, combining code from 3Dvlt_K0cFind.c and 2Dvpt_macro.c */
/* AUTHOR: Andrew Forrester <aforrester@ucla.edu> */
/* DESCRIPTION:
   * This C program returns K0c for a system of a film of liquid helium-4 (4He)...
   * The program uses Gary A. Williams's vortex-pair theory (Kosterlitz-Thouless theory) of the superfluid phase transition (i.e., its differential recursion relations) and a 4th-order Runge-Kutta integration method...

   * In the basename "2Dvpt_K0cFind", "vpt" refers to Williams's Vortex Pair Theory (of the superfluid phase transition) and "K0cFind" refers to what the program does -- it FINDs K0c.

   Inputs:
   * Kthresh     -
   * uncertainty - (the exponent must be less than -15... see "%1.15g" below... if you want all the displayed digits to be accurate)
   * DK01        - "pre-initial" and maximum adjustment of K0.  It's a "pre-initial" adjustment because we change it (divide by 2) before applying it to adjust K0.  It's the maximum adjustment of K0 because the calculated K0c (which we hope to equal the actual K0c) will be within the range K01+-DK01.  The guess for K0c (K01) must be within (actual K0c)+-DK01 to yield the correct calculate K0c.  If the calculated results for K0c are at the boundaries (either K01-DK01 or K01+DK01) then DK01 should be made larger, to find approximately where K01 should be set.  (DK01 should also be smaller than K01 so K0 is not made to go negative.)

   Output:
   * A file called "vlt_CcK01OutputDetail.dat" explains the results in detail.
   * A file called "vlt_CcK01Outputn.dat" where "n" matches the number written below.  E.g., vlt_CcK01Output03.dat.  This file writes the results in a form that can be quickly put into ... and used.
*/
/* EXT FILES: vlt_CcK01Inputn.dat */
/* COMPILE NOTES:
   * To compile, first prepare a new input file, then change the input/output filenames below,
     and then type "g++ -lm 2Dvpt_K0cFind.c" without the quotes; then to run, type "./a.out".
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
const int    N           = 2;
const double dl          = 1.0e-4;
const double Kthresh     = 5.0;
const double uncertainty = 1e-17;
const double K01         = 5;      // Initial K0, in trying to find K0c
const double DK01        = 5;      // The calculated K0c is with the range K01+-DK01



//  Function Prototypes  //
//=========================================================================//

// rk4((*f)(),*x,y[],h,n);
// EqRecRel(x,z,dzdx,n);
void rk4(void (*f)(double, double*, double*, unsigned int), double *x, double y[], double h, unsigned n);
void EqRecRel(double x, double z[], double dzdx[], unsigned n);



//  Function Definitions  //
//=========================================================================//

main(){
	int    i=0;
	char   trash[4];
	double x, l;
	double z[N], K, y;
	double Cc,K01;
	double K0,DK0,oldK,oldy;
	FILE   *outfile;

	USE 2Dvpt_K0cFind_Plot.c FOR NOW, AND IF YOU HAVE THE TIME, IMPLEMENT THE SUGGESTED IDEA(S) FOR BETTER CODE.

	// Open/name input and output files, and print headings for data to screen and output file //
	outfile = fopen("2Dvpt_K0cFind_Output.dat","w");
	fprintf(outfile,"# Source: 2Dvpt_K0cFind.c\n");
	fprintf(outfile,"# Source version: %s\n\n", "1 (2011 Dec 14 - ...)");

	// Find K0c //
	// Increment data set count (for output file) //
	i++;

	// Set K0 to initial guess K01 for K0c (which should be close to actual K0c) //
	K0 = K01;

	// Set DK0 to pre-initial value //
	DK0 = DK01;

	// Adjust K0 by finer and finer amounts until the adjustment DK0 reaches the desired precision //
	while(DK0>uncertainty){
		// Increase precision of adjustment of K0 //
		DK0 = DK0/2.0;

		// Initialize lengthscale (and energy) values //
		x = l = 0.0;
		// estimated.thrmv[2] = 0.0;

		// Initialize would-be temperature-dependent values at initial lengthscale //
		z[0] = K = K0;
		z[1] = y = exp(-PISQ*K0);
		oldK = K;
		//oldy = estimated.thrmv[1];

		// Find whether K0 should be increased or decreased by DK0 to approach K0c //
		while(1==1){
			// Step out in length scale l, calculating next K and y //
			rk4(EqRecRel, &x, z, dl, N);  // equil: rk4 2D (K,y), EqRecRel
			l += dl;
			K = z[0];
			y = z[1];

			// K is diverging, K0 too high //
			if(K>Kthresh) {
				//printf("K exploded\n");
				K0 = K0-DK0;
				break;
			}
			// K is decreasing after l=6.0, K0 too low //
			if(l>6.0 && K<oldK) {
				//printf("K decreasing\n");
				K0 = K0+DK0;
				break;
			}

			// If the above conditions have not yet been found, prepare for comparison with the next estimated state at a larger lengthscale //
			oldK = K;
			//oldy = estimated.thrmv[1];
		}

		// Bug check: Show convergence towards K0c along with approximate error DK0 //
		//printf("%1.15e %1.15e\n",K0,DK0);
	}

	// Print results to screen and output file //
	printf("%1.15g\n",K0);
	fprintf(outfile,"calculated K0c is within K01+-DK01 = %g+-%g = (%g, %g);\tcalculated K0c = %1.15g;\n", K01,DK01,K01-DK01,K01+DK01,K0);
	fflush(outfile);

	// Close output file //
	fclose(outfile);
}



// Runge-Kutta 4th order method //
void rk4(void (*f)(double, double*, double*, unsigned), double *x, double y[], double h, unsigned n){
	int i,j;
	double *k[4],*s,temp;

	for(i=0;i<4;i++)  k[i] = (double *)calloc(n,sizeof(double));
	s = (double *)calloc(n,sizeof(double));
	f(*x,y,k[0],n);
	for(i=1;i<4;i++){
		temp = h*((i+1)/2)/2;
		// using integer division trick (times double): ((i+1)/2)/2 = 0.50, 0.50, 1.00 instead of 0.75, 1.00, 1.25
		for(j=0;j<n;j++)  s[j] = y[j]+k[i-1][j]*temp;
		f(*x+h*(i/2)/2,s,k[i],n);
		// using integer division trick (times double): (i/2)/2 = 0.00, 0.50, 0.50 instead of 0.50, 0.75, 1.00
	}
	for(j=0;j<n;j++)  y[j] += h*(k[0][j]+2.0*(k[1][j]+k[2][j])+k[3][j])/6.0;
	*x += h;
	for(i=0;i<4;i++)  free((char *) k[i]);
	free((char *)s);
}



// equilibrium (K,y) recursion relations //
void EqRecRel(double x, double z[2], double dzdx[2], unsigned n){
	dzdx[0] = -B*z[0]*z[0]*z[1];
	dzdx[1] = z[1]*(4.0-2.0*PI*z[0]);
}



//  Program Notes  //
//=========================================================================//
/*

See further documentation (UCLAresearchNotes.pdf, UCLAresearchPrograms.pdf, vlt_ThermStates.c, ring1.c) for elaboration.

Dl is the increment in l.
n is the number of thermodynamic variables (n=3 given thermv[i], i=1,2,3).

DthrmvDl[i].arr[0] = dK/dl
DthrmvDl[i].arr[1] = dy/dl
DthrmvDl[i].arr[2] = de/dl

*/
