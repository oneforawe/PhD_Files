//  File Comments  //
//=========================================================================//

/* FILENAME: 3Dvlt_K0cFind2.c */
/* VERSION: 1 (2012 May 19 - ...)
            Now using a method that is more informed by the autonomous differential equations analysis, using dy/dK. */
/* AUTHOR: Andrew Forrester <aforrester@ucla.edu> */
/* DESCRIPTION:
   * ... (look in 3Dvlt_K0cFind.c to mimick notes.)

   Inputs:
   * A file called "3Dvlt_K0cFind2_Inputn.dat" where "n" matches the number written below.  E.g., vlt_CcInput3.dat.
     In the file there should be a list of values for Cc.  Cc can be any positive number, presumably.
   * ... (look in 3Dvlt_K0cFind.c to mimick notes.)

   Output:
   * ... (look in 3Dvlt_K0cFind.c to mimick notes.)
*/
/* EXT FILES: vlt_CcInputn.dat */
/* COMPILE NOTES:
   * To compile, first prepare a new input file, then change the input/output filenames below,
     and then type "g++ -lm vlt_K0cFind2.c" without the quotes; then to run, type "./a.out".
*/
/* PROGRAM IDEAS:
   * To make this code more efficient, realize that there is really only one trajectory that determines all of the K0c.  The list of Cc's can be loaded all at once, arranged into two lists (one for the positive-K-direction trajectory portion, and the other for the negative-K-direction trajectory portion), and each K0c can be calculated along the way as the trajectory is traversed.
*/

/* RESULTS:
*/



//  Function Preparation  //
//=========================================================================//

// Standard routine header files //
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Constants and labels //
const double PI     = 3.14159265358979323846;
const double PISQ   = 9.86960440108935861883;
const double PICU   = 31.00627668029981620634;

// Program parameters/inputs //
const double B      = 4.0*PICU/3.0;
const double THETA  = 0.6;
const int    N      = 1;
const double Kstar  = 0.387508189712343;
const double ystar  = 1.0/(B*Kstar);
const double c      = ystar * PISQ * ( 1.0-THETA*(1.0+log(Kstar)) );
const double Z      = Kstar * PISQ * ( 1.0-THETA*(1.0+log(Kstar)) );
const double lambda = 0.5*( 1.0 - sqrt(1.0+4.0*Z) );

const double initialAbsdK    = 1.0e-16;
const double trajectoryAbsdK = 1.0e-6;
const double uncertainty     = 1.0e-17;  // could loop, decreasing dK until within K0c is uncertainty



//  Function Prototypes  //
//=========================================================================//

// rk4((*f)(),*x,y[],h,n);
// EqRecRel(x,z,dzdx,n);
void rk4(void (*f)(double, double*, double*, unsigned int), double *x, double y[], double h, unsigned n);
void EqRecRel(double x, double z[], double dzdx[], unsigned n);



//  Function Definitions  //
//=========================================================================//

main(){
	int    nextK0cfound;
	int    count_steps=0, count_loops=0;
	char   trash[4];
	double Cc, K0c=0.0, oldK0c=-100.0;
	double dK, deltaK;
	double x, l;
	double z[N], K, y;
	double oldK, oldy;
	FILE   *infile,*outfile;

	// Open/name input and output files, and print headings for data to screen and output file //
	infile  = fopen("3Dvlt_K0cFind2_Input2.dat","r");
	outfile = fopen("3Dvlt_K0cFind2_Output2.dat","w");
	fprintf(outfile,"# Source: 3Dvlt_K0cFind2.c\n");
	fprintf(outfile,"# Source version: %s\n", "1 (2012 May 19 - ...)");
	fprintf(outfile,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g\n", PI,PISQ,PICU);
	fprintf(outfile,"# Parameter values: B=%21.21g, THETA=%g, N=%i, Kstar=%21.21g, ystar=%21.21g, c=%21.21g, Z=%21.21g, lambda=%21.21g, initialAbsdK=%21.21g, trajectoryAbsdK=%21.21g, uncertainty=%g\n", B, THETA, N, Kstar, ystar, c, Z, lambda, initialAbsdK, trajectoryAbsdK, uncertainty);
	fprintf(outfile,"# %s\t%s\n","Cc","K0c");
	printf(           "%s\t%s\n","Cc","K0c");

	// Skip first line in input file to ready data retrieval //
	fscanf(infile,"%s",trash);  // Pass "Cc"

	// Find K0c for each (Cc,K01) pair in the input file //
	while(!feof(infile)){
		// Acquire value for Cc from input file //
		fscanf(infile,"%lf",&Cc);

		// Initialize position for trajectory //
		K = Kstar;
		y = ystar;

		// Find which direction (+/-) to travel in K-direction //
		if( y < exp(-PISQ*Cc*K) )
			dK = +initialAbsdK;
		if( y > exp(-PISQ*Cc*K) )
			dK = -initialAbsdK;
		if( y== exp(-PISQ*Cc*K) ){
			K0c = Kstar;
			fprintf(outfile,"%1.15g\t%1.15g\n",Cc,K0c);
			printf(         "%1.15g\t%1.15g\n",Cc,K0c);
			continue;
		}

		// Take first step in trajectory (if not already done finding K0c) //
		oldK = K;
		oldy = y;
		K += dK;
		y += -(c/lambda)*dK;
		z[0] = y;

		// Prepare for potential additional steps in trajectory //
		if(dK<0.0) dK = -trajectoryAbsdK;
		if(dK>0.0) dK =  trajectoryAbsdK;

		// Until K0c is found within the desired precision, ... //
		while( abs(K0c-oldK0c) > uncertainty ){
			// ...cycle through checking for interception (finding K0c) and taking more steps in trajectory //
			nextK0cfound = false;
			if(dK>0.0){
				while(nextK0cfound==false){
					if( y > exp(-PISQ*Cc*K) ){
						oldK0c = K0c;
						deltaK = (exp(-PISQ*Cc*oldK)-oldy)/((y-oldy)/dK+PISQ*Cc*exp(-PISQ*Cc*oldK));
						K0c = oldK + deltaK;
						nextK0cfound = true;
					}
					else{
						oldK = K;
						oldy = y;
						rk4(EqRecRel, &K, z, dK, N);  // equil: rk4 1D (y), EqRecRel
						y = z[0];
						count_steps++;
					}
				}
			}
			if(dK<0.0){
				while(nextK0cfound==false){
					if( y < exp(-PISQ*Cc*K) ){
						oldK0c = K0c;
						deltaK = (exp(-PISQ*Cc*K)-y)/((y-oldy)/dK+PISQ*Cc*exp(-PISQ*Cc*K));
						K0c = K + deltaK;
						nextK0cfound = true;
					}
					else{
						oldK = K;
						oldy = y;
						rk4(EqRecRel, &K, z, dK, N);  // equil: rk4 1D (y), EqRecRel
						y = z[0];
						count_steps++;
					}
				}
			}
			// (in case another pass is required) //
			dK = dK/2.0;
			K = Kstar+dK;
			y = ystar-(c/lambda)*dK;
			count_loops++;
		}

		// Print results to screen and output file //
//		fprintf(outfile,"%1.15g\t%1.15g\n",Cc,K0c);
//		printf(         "%1.15g\t%1.15g\n",Cc,K0c);
		fprintf(outfile,"%1.15g\t%1.15g\n",Cc,K0c);
		printf(         "%1.15g\t%1.15g\t(steps=%i, loops=%i, dK=%g)\n",Cc,K0c,count_steps,count_loops,dK);
		fflush(outfile);
		 // Why not this instead?
		 // fprintf(outfile1,"%16.16g\t%16.16g\n",Cc,K0c);
		 // printf(          "%16.16g\t%16.16g\n",Cc,K0c);

		// Get ready for next input and search //
		K0c = 0.0;
		oldK0c = -100.0;
		count_steps = count_loops = 0;
	}

	// Close input and output files //
	fclose(infile);
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



// Equilibrium (K,y,e) recursion relations //
void EqRecRel(double x, double z[], double dzdx[], unsigned n){
	dzdx[0] = z[0]*( 6.0 - PISQ*x*(1.0-THETA*log(x)) ) / (x-B*x*x*z[0]);
}



//  Program Notes  //
//=========================================================================//
/*

See further documentation (UCLAresearchNotes.pdf, UCLAresearchPrograms.pdf, vlt_ThermStates.c, ring1.c) for elaboration.

*/
