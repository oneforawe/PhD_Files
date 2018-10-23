//  File Comments  //
//=========================================================================//

/* FILENAME: 3Dvlt_K0cFind_Plot.c */
/* VERSION: 2 (2012 Jul 27 - ...)
            B renamed as A3.
            (With the y-equation if-then statement commented out for recalculating thesis data.)
            Back to new THETA (0.6918) to see how to find a good Kthresh by plotting.
            Back to old THETA (0.6) and old Kthresh (5).
            With new recursion relations (new B and new y definition with 9 instead of 6) and new THETA (0.6918) (and trying Kthresh=10 instead of 5).
            Not allowing ac' to exceed a (See *).  (Results nearly all same as before. Cc = 0.1 gives different K0c.  See 3Dvlt_K0cFind_Output14.dat vs vlt_CcK0cOutput10.dat)
            *Everywhere with (1-THETA*log(K)) set log(K) to zero if it goes positive.
            C is now out of the recursion relations, normal settings */
/* AUTHOR: Andrew Forrester <aforrester@ucla.edu> */
/* DESCRIPTION:
   * This C program takes (one or more values of) Cc and returns the corresponding value(s) of K0c for a system of liquid helium-4 (4He).  More description of Cc, K0c, and the system follows.
   * The system is a spherical region (of size greater than the smallest relevant length scale a0, the diameter of the smallest vortex loops, where the size a0 itself, in terms of meters, say, is irrelevant) of liquid 4He at some constant pressure, near or below the critical transition temperature Tc (which depends on pressure) separating normal fluid and superfluid (aka the lambda temperature).
   * The program calculates K0c, the critical "bare" superfluid ratio at the smallest relevant length scale a0, for the system, given that Cc, the critical dimensionless vortex-loop-core energy, takes on some specified value(s).  The value of Cc also happens to specify the pressure of the system.
   * The program uses Gary A. Williams's vortex-loop theory of the superfluid phase transition (i.e., its differential recursion relations) and a 4th-order Runge-Kutta integration method...

   * In the basename "vlt_K0cFind", "vlt" refers to Williams's Vortex Loop Theory (of the superfluid phase transition) and "K0cFind" refers to what the program does -- it FINDs K0c.

   Inputs:
   * A file called "vlt_CcK01Inputn.dat" where "n" matches the number written below.  E.g., vlt_CcK01Input3.dat.
     In the file there should be a list of 2-tuples (Cc and K01).  Cc can be anything, presumably.  The corresponding K01 value must be close to the actual corresponding K0c for the program to find the actual K0c.  If you don't know the approximate value of the actual K0c, you can try a number and see whether the program "max's out" (meaning your guess was too low) or "min's out" (meaning your guess was too high).
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
     and then type "g++ -lm vlt_K0cFind.c" without the quotes; then to run, type "./a.out".
*/

/* RESULTS:
       Given  Cc = 1.03  and  K01 = 0.3,
       Previous calculated values have been K0c = 0.309146998265557
                                                  0.3091469984282580
                                                  0.3091469984216758
       The most updated calculation gets    K0c = 0.3091469984282583
        and now the value (the last few digits) of K0c should not
        depend on K01, since the parameter "uncertainty" is low enough.

       See files such as vlt_CcK0cOutput01.dat for further results.


Further Results:
Cc	K0c(Cc)
1.20	0.2805443089652032
1.10	0.2964689389987730
1.06	0.3035407129346564
1.05	0.3053794236011743
1.04	0.3072479227697506
1.03	0.3091469984216724
1.02	0.3110774675059972
1.01	0.3130401773072950
1.00	0.3150360068926113
0.99	0.3170658686428607
0.98	0.3191307098748213
0.97	0.3212315145597378
0.90	0.3370353710663317
0.80	0.3636228655823348
0.70	0.3965669728967915
0.60	0.4387310072100324
0.55	0.4646970470439853
0.50	0.4951082899772941
0.40	0.5753434296585114
0.30	0.7010909420810506
0.20	0.9346549781449002
0.10	1.576187055119773
0.00	5.105045166053435

(See RhosAmps.ods for Cc(P), derived from Ahlers' 1973 data):
P (bar)	Cc(P)	K0c(Cc) (from K0cFind.c)
um...
0.050	1.105023144261770	0.295611433115442
1.646	1.058195551375160	0.303870330660806
7.328	0.905924430200288	0.335617836037364
15.031	0.730725772346878	0.385625265795951
18.180	0.667963030000677	0.408889237394402
22.533	0.588565898089815	0.444331744979223
25.868	0.532988906023542	0.474486704164683
wait, isn't it this:
0.050	1.10500797392886
1.646	1.0581889918361
7.328	0.905918431108949
15.031	0.730685115724445
18.180	0.667831302441093
22.533	0.588281338805611
25.868	0.532635362503492
'cause this is in 3Dvlt_quench.c:
	P =  0.050;	Cc = 1.10500797392886;	K0c = 0.295614012972922;
*/



//  Function Preparation  //
//=========================================================================//

// Standard routine header files //
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Constants and labels //
const double PI    = 3.14159265358979323846;
const double PISQ  = 9.86960440108935861883;
const double PICU  = 31.00627668029981620634;

// Program parameters/inputs //
const int    OpA  = 1;   //  Option A   Set recursion relation behavior
                         //              1 -> "itlog flow": use if-then statement to limit ac at a -- "if-then log flow" of EqRecRel1;
                         //              2 -> "w/log flow": no if-then, using logarithm so ac may exceed a (two-fixed points) -- "with log flow" of EqRecRel2;
const double A3    = 4.0*PICU/3.0;
//const double A3    = 2.0*PICU*PISQ/3.0;
const double THETA = 0.6;
//const double THETA = 0.6918;
const int    N     = 2;
const double dl           = 1.0e-4;
const int    LDataPerOne  = 30;      // Number of l data points collected per Delta l = 1 (from l=0 to l=1, or l=1 to l=2, etc.)
const int    MRecDataSets = 10;      // Maximum number of recorded data sets
const double uncertainty  = 1.0e-17;
const double DK01         = 2.0;     // The calculated K0c is with the range K01+-DK01



//  Function Prototypes  //
//=========================================================================//

// rk4((*f)(),*x,y[],h,n);
// EqRecRel(x,z,dzdx,n);
void rk4(void (*f)(double, double*, double*, unsigned int), double *x, double y[], double h, unsigned n);
void EqRecRel1(double x, double z[], double dzdx[], unsigned n);
void EqRecRel2(double x, double z[], double dzdx[], unsigned n);
void PressEnterToContinue();



//  Function Definitions  //
//=========================================================================//

int main(){
	int    i,j,k,m;
	char   trash[4];
	double x, l;
	double z[N], K, y;
	double oldK, oldy;
	double slope;
	double Cc,K01;
	double K0,DK0;
	int    LSkip = 1.0/dl/LDataPerOne;
	int    DataSets = ceil(log(DK01/uncertainty)/log(2.0));
	int    DataSkip = ceil((DataSets*1.0)/(MRecDataSets*1.0));
	int    RecDataSets = floor(DataSets/DataSkip);
	FILE   *infile,*outfile1,*outfile2,*outfile3;
	char   *outfilename1,*outfilename2,*outfilename3;
	void   (*EqRecRel)(double x, double z[], double dzdx[], unsigned n);

	switch(OpA){ // Set recursion relation behavior
		case 1:  // "itlog flow": use if-then statement to limit ac at a -- "if-then log flow"
			 EqRecRel = EqRecRel1;
			 break;
		case 2:  // "w/log flow": no if-then, using logarithm so ac may exceed a (two-fixed points) -- "with log flow"
			 EqRecRel = EqRecRel2;
			 break;
	}

	/* Warning to check how many data sets will be recorded */
	printf("\nThere will be %i data sets, and %i of them will be recorded.\n", DataSets, RecDataSets );
	printf("If that sounds like too many,\n");
	printf("Press Ctrl+C now to cancel,\n");
	printf("otherwise\n");
	PressEnterToContinue();

	// Open/name input and output files, and print headings for data to screen and output file //
	infile   = fopen("3Dvlt_K0cFind_Plot_Input5.dat","r");               // change number here  (when using new Input file)
	asprintf(&outfilename1, "3Dvlt_K0cFind_Plot_Output5_Op%i.dat", OpA); // and here, and below (for outfile3)
	asprintf(&outfilename2, "3Dvlt_K0cFind_Plot_OutputDetail.dat", OpA);
	outfile1 = fopen(outfilename1,"w");  //  E.g., "3Dvlt_K0cFind_Plot.out"
	outfile2 = fopen(outfilename2,"w");  //  E.g., "3Dvlt_K0cFind_Plot.out"
	fprintf(outfile1,"# Filename: %s\n", outfilename1);
	fprintf(outfile1,"# Source: %s\n", "3Dvlt_K0cFind_Plot.c");
	fprintf(outfile1,"# Source version: %s\n", "2 (2012 Jul 27 - ...)");
	fprintf(outfile1,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g\n", PI,PISQ,PICU);
	fprintf(outfile1,"# Parameter values: A3=%21.21g, THETA=%g, N=%i, dl=%g, LDataPerOne=%i, MRecDataSets=%i, uncertainty=%g, DK01=%g\n", A3, THETA, N, dl, LDataPerOne, MRecDataSets, uncertainty, DK01);
	fprintf(outfile2,"# Filename: %s\n", outfilename2);
	fprintf(outfile2,"# Source: %s\n", "3Dvlt_K0cFind_Plot.c");
	fprintf(outfile2,"# Source version: %s\n", "2 (2012 Jul 27 - ...)");
	fprintf(outfile2,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g\n", PI,PISQ,PICU);
	fprintf(outfile2,"# Parameter values: A3=%21.21g, THETA=%g, N=%i, dl=%g, LDataPerOne=%i, MRecDataSets=%i, uncertainty=%g, DK01=%g\n", A3, THETA, N, dl, LDataPerOne, MRecDataSets, uncertainty, DK01);
	// fprintf(outfile,"%s\t%s\n","Cc","K0c");
	printf(         "%s\t%s\n","Cc","K0c");

	// Skip first line in input file to ready data retrieval //
	fscanf(infile,"%s",trash);  // Pass "Cc"
	fscanf(infile,"%s",trash);  // Pass "K01"

	// Find K0c for each (Cc,K01) pair in the input file //
	i = 0;
	while(!feof(infile)){

		// Acquire values for Cc and K01 from input file //
		fscanf(infile,"%lf%lf",&Cc,&K01);
		if(feof(infile)) break;
		// Increment Cc,K0c data-set/case count (for output file) //
		i++;
		m=1;  // reset l-data set counter for current Cc,K0c data set

		// Prepare flow-plotting output file //
		//asprintf(&filename, "3Dvlt_K0cFind_Plot_Cc%g_K01%g.dat", Cc,K01);
		asprintf(&outfilename3, "3Dvlt_K0cFind_Plot_Output5_flow_Cc%3.2f_K01%g.dat", Cc,K01); // change number here
		outfile3 = fopen(outfilename3,"w");  //  E.g., "3Dvlt_K0cFind_Plot.out"
		fprintf(outfile3,"# Filename: %s\n", outfilename1);
		fprintf(outfile3,"# Source: %s\n", "3Dvlt_K0cFind_Plot.c");
		fprintf(outfile3,"# Source version: %s\n", "2 (2012 Jul 27 - ...)");
		fprintf(outfile3,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g\n", PI,PISQ,PICU);
		fprintf(outfile3,"# Parameter values: A3=%21.21g, THETA=%g, N=%i, dl=%g, LDataPerOne=%i, MRecDataSets=%i, uncertainty=%g, DK01=%g\n", A3, THETA, N, dl, LDataPerOne, MRecDataSets, uncertainty, DK01);
		fprintf(outfile3,"# %s\t%s\t%s\n\n", "l", "K", "y");

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
			z[1] = y = exp(-PISQ*K0*Cc);
			if(j%DataSkip==0 || DK0/2.0<uncertainty){
				fprintf(outfile3,"# DK0 = %g\t(data set %i)\n", DK0,m);  m++;
				fprintf(outfile3,"%g\t%17.17g\t%17.17g\n", l, K, y);
			}

			// Find whether K0 should be increased or decreased by DK0 to approach K0c //
			k = 0;
			slope = 1.0;
			while(slope>0.0){
				// Step out in length scale l, calculating next K and y //
				oldK = K;
				oldy = y;
				rk4(EqRecRel, &x, z, dl, N);  // equil: rk4 2D (K,y), EqRecRel
				l += dl;
				k++;       // step in l
				K = z[0];
				y = z[1];
				slope = (y-oldy)/(K-oldK);

				// Print data points to output file //
				if( (j%DataSkip==0 || DK0/2.0<uncertainty) && k%LSkip==0)
					fprintf(outfile3,"%g\t%17.17g\t%17.17g\n", l, K, y);
			}
			if(K<oldK){ // then K0 was too low
				K0 = K0+DK0;
				if(y<oldy) { printf("Error: Expected the trajectory to up-left (when K0 is too low and after swooping by the fixed point). This was not the case."); exit(EXIT_FAILURE); }
			}
			if(K>oldK){ // then K0 was too high
				K0 = K0-DK0;
				if(y>oldy) { printf("Error: Expected the trajectory to down-right (when K0 is too high and after swooping by the fixed point). This was not the case."); exit(EXIT_FAILURE); }
			}
			if(K==oldK){
				if(y<oldy){ // then K0 was too high
					K0 = K0-DK0;
				}
				if(y>oldy){ // then K0 was too low
					K0 = K0+DK0;
				}
			}

			// Comment... //
			if(j%DataSkip==0 || DK0/2.0<uncertainty)  fprintf(outfile3,"\n");
			if(DK0/2.0<uncertainty)  break;  // Is this necessary?

			// Show convergence towards K0c along with approximate error DK0 (for bug checking) //
			printf("step %i: %1.17e %1.17e\n", j,K0,DK0);
		}

		// Print results to screen and output file //
		//printf(          "%1.17g\t%1.17g\n", Cc,K0);
		printf(          "%3.2f\t%1.17g\n", Cc,K0);
		//fprintf(outfile1,"%1.17g\t%1.17g\n", Cc,K0);
		//fprintf(outfile1,"case %i:\t\tCc = %1.15g;\tK0c = %1.15g;\tbreak;\n", i,Cc,K0);
		//fprintf(outfile1,"case %i:\t\tCc = %1.15g;\tK0c = %17.17g;\tbreak;\n", i,Cc,K0);
		fprintf(outfile1,"\t\t\tcase %i:\t\tCc = %3.2f;\tK0c = %17.17g;\tbreak;  // ac/a = K^theta*exp(C) = %3f (@ Tc, l=0)\n", i,Cc,K0, pow(K0,THETA)*exp(Cc) );
		fflush(outfile1);
		fprintf(outfile2,"case %i:\tCc = %1.15g;\tcalculated K0c is within K01+-DK01 = %g+-%g = (%g, %g);\tcalculated K0c = %1.15g;\n", i,Cc,K01,DK01,K01-DK01,K01+DK01,K0);
		fflush(outfile2);
		
		// Close flow-plotting output file for certain Cc,K01 pair //
		fclose(outfile3);
	}

	// Close input and output files //
	fclose(infile);
	fclose(outfile1);
	fclose(outfile2);
	
	return EXIT_SUCCESS;
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



// Equilibrium (K,y) recursion relations (version 1, proper) //
void EqRecRel1(double x, double z[], double dzdx[], unsigned n){
	dzdx[0] = z[0]-A3*z[0]*z[0]*z[1];
	if( log(z[0]) > 0.0 ) // if ac' > a, set ac' = a (replace log(K[i]) with zero)
		dzdx[1] = z[1]*( 6.0 - PISQ*z[0] );
	else
		dzdx[1] = z[1]*( 6.0 - PISQ*z[0]*(1.0-THETA*log(z[0])) );
//	dzdx[2] = -PI*exp(-3.0*x)*z[1];
}



// Equilibrium (K,y) recursion relations (version 2) //
void EqRecRel2(double x, double z[], double dzdx[], unsigned n){
	dzdx[0] = z[0]-A3*z[0]*z[0]*z[1];
	//if( log(z[0]) > 0.0 ) // if ac' > a, set ac' = a (replace log(K[i]) with zero)
	//	dzdx[1] = z[1]*( 6.0 - PISQ*z[0] );
	//else
		dzdx[1] = z[1]*( 6.0 - PISQ*z[0]*(1.0-THETA*log(z[0])) );
//	dzdx[2] = -PI*exp(-3.0*x)*z[1];
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
