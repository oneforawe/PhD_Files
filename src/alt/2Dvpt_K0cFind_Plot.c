//  File Comments  //
//=========================================================================//

/* FILENAME: 2Dvpt_K0cFind_Plot.c */
/* VERSION: 1 (2011 Dec 14 - ...)
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
/* NOTE IDEA(S) FOR BETTER CODE:
   * For the "If y is decreasing" block of code, test instead whether the logarithm of y has the same slope over Delta l of about 1000.
   When y is actually settled into a decrease to zero (and not preparing eventually to swing back up and explode), y takes on an exponential form whose slope remains the same y - y0*exp[(4-2pi*K)l]
   * Or, EVEN BETTER, if K stops approaching 2.0/pi (i.e., perhaps determined when decrease in K, K{i}-K{i-1}, is less than the difference between K and K*, K{i}-2.0/pi), then K0 was too high, and if K goes below 2.0/pi, the K0 was too low.  (With this test "lcutoff" will determine when the loop stops rather than "uncertainty".)
*/
/* SOME RESULTS:

   First, I was using  K0c = 0.747853  from Gary's calculations.
   lcutoff = 4000.0  yields about  0.7478523616
   lcutoff = 4500.0  yields  K0c = 0.74785238
   lcutoff = 4700.0  yields  K0c = 0.74785238
   lcutoff = 4750.0  yields  K0c = 0.74785238
 
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
const int    N            = 2;
const double dl           = 1.0e-4;
const double lcutoff      = 4750.0;  // If y hasn't exploded yet, test here whether y is increasing or decreasing (must be quite large to get K0c precisely)
const int    LDataPerOne  = 1;     // Number of l data points collected per Delta l = 1 (from l=0 to l=1, or l=1 to l=2, etc.)
const int    MRecDataSets = 10;      // Maximum number of recorded data sets
//const double Kthresh      = 3.0;
const double ythresh      = 1.0;
const double uncertainty  = 1e-11;    //  1e-17
//const double K01          = 0.74785;      // Initial K0, in trying to find K0c
//const double DK01         = 0.74785;      // The calculated K0c is with the range K01+-DK01
const double K01          = 0.5;      // Initial K0, in trying to find K0c
const double DK01         = 0.5;      // The calculated K0c is with the range K01+-DK01



//  Function Prototypes  //
//=========================================================================//

// rk4((*f)(),*x,y[],h,n);
// EqRecRel(x,z,dzdx,n);
void rk4(void (*f)(double, double*, double*, unsigned int), double *x, double y[], double h, unsigned n);
void EqRecRel(double x, double z[], double dzdx[], unsigned n);



//  Function Definitions  //
//=========================================================================//

main(){
	int    i=0,j,k,m=1,wentup=false;
	char   trash[4];
	double x, l;
	double z[N], K, y;
	double K0,DK0,oldK,oldy;
	int    LSkip = 1.0/dl/LDataPerOne;
	int    DataSets = ceil(log(DK01/uncertainty)/log(2.0));
	int    DataSkip = ceil((DataSets*1.0)/(MRecDataSets*1.0));
	int    RecDataSets = floor(DataSets/DataSkip);
	FILE   *outfile,*outfile2;
	char   *filename,*filename2;

	/* Warning to check how many data sets will be recorded //
	printf("\nThere will be %i data sets, and %i of them will be recorded.", DataSets, RecDataSets );
	printf("\nIf that sounds like too many, hit Ctrl+C now, otherwise...\n");
	PressEnterToContinue();*/

	// Open/name input and output files, and print headings for data to screen and output file //
	asprintf(&filename,  "2Dvpt_K0cFind_Output2.dat");
	asprintf(&filename2, "2Dvpt_K0cFind_Plot_Output.dat");
	outfile  = fopen(filename, "w");
	outfile2 = fopen(filename2,"w");
	fprintf(outfile, "# Filename: %s\n", filename);
	fprintf(outfile2,"# Filename: %s\n", filename2);
	fprintf(outfile, "# Source: 2Dvpt_K0cFind_Plot.c\n");
	fprintf(outfile2,"# Source: 2Dvpt_K0cFind_Plot.c\n");
	fprintf(outfile, "# Source version: %s\n", "1 (2011 Dec 14 - ...)");
	fprintf(outfile2,"# Source version: %s\n", "1 (2011 Dec 14 - ...)");
	fprintf(outfile2,"# %s\t%s\t%s\n\n", "l", "K", "y");
	printf("K0\ty0\n");

	// Find K0c //
	// Set K0 to initial guess K01 for K0c (which should be close to actual K0c) //
	K0 = K01;

	// Set DK0 to pre-initial value //
	DK0 = DK01;

	// Adjust K0 by finer and finer amounts until the adjustment DK0 reaches the desired precision //
	j=-1;
	while(DK0>uncertainty){
		// Increase precision of adjustment of K0 //
		j++; // step in precision
		DK0 = DK0/2.0;

		// Initialize lengthscale (and energy) values //
		x = l = 0.0;
		// estimated.thrmv[2] = 0.0;

		// Initialize would-be temperature-dependent values at initial lengthscale //
		z[0] = K = K0;
		z[1] = y = exp(-PISQ*K0);
		oldy = y;
		//oldK = K;
		printf("%10.10g\t%10.10g\t",z[0],z[1]);
		if(K0==K01)
			printf("K0 = K01\n");
		else{
			if(wentup==true)
				printf("K0 went up\n");
			else
				printf("K0 went down\n");
		}

		// Find whether K0 should be increased or decreased by DK0 to approach K0c //
		if(j%DataSkip==0 || DK0/2<uncertainty){
			fprintf(outfile2,"# DK0 = %g\t(data set %i)\t", DK0,m);  m++;
			if(K0==K01)
				fprintf(outfile2,"K0 = K01\n");
			else{
				if(wentup==true)
					fprintf(outfile2,"K0 went up\n");
				else
					fprintf(outfile2,"K0 went down\n");
			}
			fprintf(outfile2,"%g\t%g\t%g\n", l, K, y);
		}
		k=0;
		while(y!=NAN){
			// Step out in length scale l, calculating next K and y //
			rk4(EqRecRel, &x, z, dl, N);  // equil: rk4 2D (K,y), EqRecRel
			l += dl;
			k++;       // step in l
			K = z[0];
			y = z[1];

			// Print data points to output file //
			if( (j%DataSkip==0 || DK0/2<uncertainty) && k%LSkip==0)  fprintf(outfile2,"%g\t%g\t%g\n", l, K, y);

			// If y is diverging, K0 too low //
			if( y>ythresh || (l>lcutoff && y>oldy) ) {
				//printf("y exploded\n");
				K0 = K0+DK0;
				wentup = true;
				break;
			}
			// If y is decreasing, K0 too high (K0c yields the slowest decrease) //
			if(l>lcutoff && y<oldy) {
				//printf("y decreasing too fast\n");
				K0 = K0-DK0;
				wentup = false;
				break;
			}

			// If the above conditions have not yet been found, prepare for comparison with the next estimated state at a larger lengthscale //
			oldy = y;
			//oldK = K;
		}
		if(j%DataSkip==0 || DK0/2<uncertainty)  fprintf(outfile2,"\n");
		if(DK0/2<uncertainty)  break;

		// Bug check: Show convergence towards K0c along with approximate error DK0 //
		//printf("%1.15e %1.15e\n",K0,DK0);
	}

	// Print results to screen and output file //
	printf("K0c = %10.10g\n",K0);
	fprintf(outfile,"calculated K0c is within K01+-DK01 = %g+-%g = (%g, %g);\tcalculated K0c = %10.10g;\n", K01,DK01,K01-DK01,K01+DK01,K0);
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



void PressEnterToContinue(){
	int c;
	printf("Press ENTER to continue... ");
	fflush(stdout);
	do c = getchar(); while ((c != '\n') && (c != EOF));
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
