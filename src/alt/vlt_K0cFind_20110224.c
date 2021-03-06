/*  File Comments  */
/*=========================================================================*/

/* This is the version right before I added Cc in the recursion relations and changed the code to make it possible to pass Cc to vltRecRel. */

/* FILENAME: vlt_K0cFind.c */
/* VERSION: 4 (2010 07 28 - 2011 02 24)
/* AUTHOR: Andrew Forrester <aforrester@ucla.edu>  */
/* DESCRIPTION:
   * This C program takes (one or more values of) Cc and returns the corresponding value(s) of K0c for a system of liquid helium-4 (4He).  More description of Cc, K0c, and the system follows.
   * The system is a spherical region (of size greater than the smallest relevant length scale a0, the diameter of the smallest vortex loops, where the size itself is irrelevant) of liquid 4He at some constant pressure, near or below the critical transition temperature Tc (which depends on pressure) separating normal fluid and superfluid (aka the lambda temperature).
   * The program calculates K0c, the critical "bare" superfluid ratio at the smallest relevant length scale a0, for the system, given that Cc, the critical dimensionless vortex-loop-core energy, takes on some specified value(s).  The value of Cc also happens to specify the pressure of the system.
   * The program uses Gary A. Williams's vortex-loop theory of the superfluid phase transition (i.e., its differential recursion relations) and a 4th-order Runge-Kutta integration method...

   * In the basename "vlt_K0cFind", "vlt" refers to Williams's Vortex Loop Theory (of the superfluid phase transition) and "K0cFind" refers to what the program does -- it FINDs K0c.

   Inputs:
   * A file called "vlt_CcK01Inputn.dat" where "n" matches the number written below.  E.g., vlt_CcK01Input3.dat.
     In the file there should be a list of 2-tuples (Cc and K01).  Cc can be anything, presumably. The corresponding K01 value must be close to the actual corresponding K0c for the program to find the actual K0c.  If you don't know the approximate value of the actual K0c, you can try a number and see whether the program "max's out" (meaning your guess was too low) or "min's out" (meaning your guess was too high).
   * Kthresh     -
   * uncertainty - (the exponent must be less than -15... see "%1.15g" below... if you want all the displayed digits to be accurate)
   * DK01        - "pre-initial" adjustment of K0.  It's a "pre-initial" adjustment because we change it (divide by 2) before using it.  This determines the initial adjustment in K0 and thus sets the requirement for how close the guess for K0c (K01) has to be to the actual K0c.  Essentially, K01 must be within K0c +- DK0.

   Output:
   * A file called "vlt_CcK01Outputn.dat" where "n" matches the number written below.  E.g., vlt_CcK01Output3.dat.
*/
/* EXT FILES: vlt_derk.c */
/* COMPILE NOTES:
   * To compile, type "gcc -lm vlt_K0cFind.c vlt_derk.c" without the quotes; then to run, type "./a.out".
   * NOTE: Be sure that N=2 in vlt_derk.c when compiling.
*/

/* RESULTS:
       Given  Cc = 1.03  and  K01 = 0.3,
       Previous calculated values have been K0c = 0.309146998265557
                                                  0.3091469984282580
                                                  0.3091469984216758
       The most updated calculation gets    K0c = 0.3091469984282583
        and now the value (the last few digits) of K0c should not
        depend on K01, since the parameter "uncertainty" is low enough.

       See files such as vlt_CcK0cOutput1.dat for further results.


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
0.050	1.105023144261770	0.295611433115442
1.646	1.058195551375160	0.303870330660806
7.328	0.905924430200288	0.335617836037364
15.031	0.730725772346878	0.385625265795951
18.180	0.667963030000677	0.408889237394402
22.533	0.588565898089815	0.444331744979223
25.868	0.532988906023542	0.474486704164683
*/



/*  Function Preparation  */
/*=========================================================================*/

/* Standard routine header files */
#include <stdio.h>
#include <math.h>

/* Constants and labels */
#define IN          // input label
#define OUT         // output label
#define INOUT       // input/output label
#define PI          3.14159265358979323846
#define PISQ        9.86960440108935861883
#define A           4.0*PI*PI*PI/3.0
#define N           2

/* Program parameters/inputs */
#define Kthresh     5.0
#define uncertainty 1e-17
#define DK01        0.1

/* Data type definitions */
typedef struct {
	double l;
	double tempv;
	double thrmv[N];
} STATE;

typedef struct {
	double arr[N];
} RETARRAY;



/*  Function Prototypes  */
/*=========================================================================*/
//   derk(func(),calculated,Dl,n)
//   ringrecrel(l,z)

extern STATE derk(RETARRAY (*func)(double l, double z[]), STATE calculated, double Dl, unsigned n);

RETARRAY vltRecRel(double l, double z[]);



/*  Function Definitions  */
/*=========================================================================*/

main(){
	int i=0;
	char trash[4];
	double Cc,K01,K0,DK0,oldK,oldy;
	STATE estimated;
	FILE *infile,*outfile;

	/* Open/name input and output files, and print headings for data to screen and output file */
	infile  = fopen("vlt_CcK01Input5.dat","r");
	outfile = fopen("vlt_CcK0cOutput5.dat","w");
	// fprintf(outfile,"%s\t%s\n","Cc","K0c");
	printf(         "%s\t%s\n","Cc","K0c");

	/* Skip first line in input file to ready data retrieval */
	fscanf(infile,"%s",trash);  // Pass "Cc"
	fscanf(infile,"%s",trash);  // Pass "K01"

	/* Find K0c for each (Cc,K01) pair in the input file */
	while(!feof(infile)){
		/* Increment data set count (for output file) */
		i++;

		/* Acquire values for Cc and K01 from input file */
		fscanf(infile,"%lf%lf",&Cc,&K01);

		/* Set K0 to initial guess K01 for K0c (which should be close to actual K0c) */
		K0 = K01;

		/* Set DK0 to pre-initial value */
		DK0 = DK01;

		/* Adjust K0 by finer and finer amounts until the adjustment DK0 reaches the desired precision */
		while(DK0>uncertainty){
			/* Increase precision of adjustment of K0 */
			DK0 = DK0/2;

			/* Initialize lengthscale (and energy) values */
			estimated.l = 0.0;
			// estimated.thrmv[2] = 0.0;

			/* Initialize would-be temperature-dependent values at initial lengthscale */
			estimated.thrmv[0] = K0;
			estimated.thrmv[1] = 1.0/(exp(PISQ*K0*Cc)-1.0);  // exp(-PISQ*K0*Cc);
			oldK = estimated.thrmv[0];
			//oldy = estimated.thrmv[1];

			/* Find whether K0 should be increased or decreased by DK0 to approach K0c */
			while(1==1){
				/* Calculate next estimated state, incrementing lengthscale */
				estimated = derk(IN vltRecRel, INOUT estimated, IN 0.0001, IN N);  // Dl=0.0001, n=N=2

				/* K is diverging, K0 too high */
				if(estimated.thrmv[0]>Kthresh) {
					//printf("K exploded\n");
					K0 = K0-DK0;
					break;
				}
				/* K is decreasing after l=6.0, K0 too low */
				if(estimated.l>6.0 && estimated.thrmv[0]<oldK) {
					//printf("K decreasing\n");
					K0 = K0+DK0;
					break;
				}

				/* If the above conditions have not yet been found, prepare for comparison with the next estimated state at a larger lengthscale */
				oldK = estimated.thrmv[0];
				//oldy = estimated.thrmv[1];
			}

			/* Show convergence towards K0c along with approximate error DK0 (for bug checking)
			printf("%1.15e %1.15e\n",K0,DK0); */
		}

		/* Print results to screen and output file */
		printf("%1.15g\t%1.15g\n",Cc,K0);
		fprintf(outfile,"%s%i%s\t\t%s%1.15g%s\t%s%1.15g%s\t%s\n","case ",i,":","Cc = ",Cc,";","K0c = ",K0,";","break;");
		//fprintf(outfile,"%1.15g\t%1.15g\n",Cc,K0);
		 // Why not this instead?
		 // printf("%16.16g\t%16.16g\n",Cc,K0);
		 // fprintf(outfile,"%16.16g\t%16.16g\n",Cc,K0);
	}

	/* Close input and output files */
	fclose(infile);
	fclose(outfile);
}



// vortex loop theory recursion relations
RETARRAY vltRecRel(double l, double z[]){
	RETARRAY dzdl;

	dzdl.arr[0] = z[0]-A*z[1]*z[0]*z[0];
	dzdl.arr[1] = z[1]*(6.0-PISQ*z[0]*(1.0-0.6*log(z[0])));
	// dzdl.arr[2] = -PI*z[1]*exp(-3.0*l);

	return dzdl;
}
