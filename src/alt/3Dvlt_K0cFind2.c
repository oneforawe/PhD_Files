//  File Comments  //
//=========================================================================//

/* FILENAME: 3Dvlt_K0cFind2.c */
/* VERSION: 2 (2012 May 28 - ...)
            Now using a GSL (GNU Science Library) subroutine to find the intersection of the exponential and the trajectory.
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

// Constants and labels //
const double PI     = 3.14159265358979323846;
const double PISQ   = 9.86960440108935861883;
const double PICU   = 31.00627668029981620634;

// Program parameters/inputs //
const double B      = 4.0*PICU/3.0;
const double THETA  = 0.6;
const int    N      = 1;
//const double Kstar  = 0.387508189712343;     //  from 3Dvlt_KstarFind.c
const double Kstar  = 0.387508189712343409;  //  Kstar "2" from 3Dvlt_KstarFind.c
//const double Kstar  = 4.146766416946759293;  //  Kstar "1" from 3Dvlt_KstarFind.c

const double ystar  = 1.0/(B*Kstar);
const double c      = ystar * PISQ * ( 1.0-THETA*(1.0+log(Kstar)) );
const double Z      = Kstar * PISQ * ( 1.0-THETA*(1.0+log(Kstar)) );
const double lambda = 0.5*( 1.0 - sqrt(1.0+4.0*Z) );

const int    Op = 1;    //  Option (1 -> const. dK step size; 2 -> const. dl step size)
const double initialAbsdK    = 1.0e-6;  // first step size, using eigenvector direction
const double trajectoryAbsdK = 1.0e-6;  //  second and later step size taken in trajectory, using scaling/recursion equations
//const double initialAbsdK    = 1.0e-16;
//const double trajectoryAbsdK = 8.0e-9;
const double uncertainty     = 1.0e-2;  // could loop, decreasing dK until within K0c is uncertainty

struct TrajectoryIntercept_params {
	double Ccval, absdK, yleft, yright, Kleft, Kright;
};



//  Function Prototypes  //
//=========================================================================//

// rk4((*f)(),*x,y[],h,n);
// EqRecRel(x,z,dzdx,n);
void rk4 (void (*f)(double, double*, double*, unsigned int), double *x, double y[], double h, unsigned n);
void EqRecRel (double x, double z[], double dzdx[], unsigned n);
double FindTrajectoryIntercept (double CcVal, double AbsdK, double yLeft, double yRight, double KLeft, double KRight);
double TrajectoryIntercept (double x, void *params);



//  Function Definitions  //
//=========================================================================//

main(){
	int    nextK0cfound;
	int    count_steps=0, count_loops=0;
	char   trash[4];
	double Cc, K0c=0.0, oldK0c=-100.0;
	double dKi, dK, deltaK;
	double x, l;
	double z[N], K, y;
	double oldK, oldy;
	FILE   *infile,*outfile;

	// Open/name input and output files, and print headings for data to screen and output file //
	infile  = fopen("3Dvlt_K0cFind2_Input3.dat","r");
	outfile = fopen("3Dvlt_K0cFind2_Output3.dat","w");
	fprintf(outfile,"# Source: 3Dvlt_K0cFind2.c\n");
	fprintf(outfile,"# Source version: %s\n", "2 (2012 May 28 - ...)");
	fprintf(outfile,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g\n", PI,PISQ,PICU);
	fprintf(outfile,"# Parameter values: B=%21.21g, THETA=%g, N=%i, Kstar=%21.21g, ystar=%21.21g, c=%21.21g, Z=%21.21g, lambda=%21.21g, dataskip=%i, Op=%i, initialAbsdK=%21.21g, trajectoryAbsdK=%21.21g, uncertainty=%g\n", B, THETA, N, Kstar, ystar, c, Z, lambda, dataskip, Op, initialAbsdK, trajectoryAbsdK, uncertainty);
	fprintf(outfile,"# %s\t%s\n","Cc","K0c");
	printf(           "%s\t%s\n","Cc","K0c");
	fflush(outfile);

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
		if( y < exp(-PISQ*Cc*K) ){
			dKi = +initialAbsdK;
			dK  = +trajectoryAbsdK;
		}
		if( y > exp(-PISQ*Cc*K) ){
			dKi = -initialAbsdK;
			dK  = -trajectoryAbsdK;
		}
		if( y== exp(-PISQ*Cc*K) ){
			K0c = Kstar;
			fprintf(outfile,"%1.15g\t%1.15g\n",Cc,K0c);
			printf(         "%1.15g\t%1.15g\n",Cc,K0c);
			fflush(outfile);
			continue;
		}

		// Take first step in trajectory (if not already done finding K0c) //
		oldK = K;
		oldy = y;
		K += dKi;
		y += -(c/lambda)*dKi;
		z[0] = y;

		// Until K0c is found within the desired precision, ... //
		while( fabs(K0c-oldK0c) > uncertainty ){
			printf(" K0c hasn't reached the right precision yet: K0c=%1.15g\t( fabs(K0c-oldK0c) = [%g] > [%g] = uncertainty;  dK=[%g] )\n", K0c, fabs(K0c-oldK0c), uncertainty, dK);

			// ...cycle through checking for interception (finding K0c) and taking more steps in trajectory //
			nextK0cfound = false;
			if(dK>0.0){
				while(nextK0cfound==false){
					if( y > exp(-PISQ*Cc*K) ){
						oldK0c = K0c;
						deltaK = FindTrajectoryIntercept (Cc, fabs(dK), oldy, y, oldK, K);
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
						deltaK = FindTrajectoryIntercept (Cc, fabs(dK), y, oldy, K, oldK);
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
			K = Kstar + dKi;
			y = ystar - (c/lambda)*dKi;
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



double
FindTrajectoryIntercept (double CcVal, double AbsdK, double yLeft, double yRight, double KLeft, double KRight) {
	int status;
	int iter = 0, max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r = 0.0, r_expected = AbsdK/2.0;
	double x_lo = 0.0, x_hi = AbsdK;
	gsl_function F;
	struct TrajectoryIntercept_params params = {CcVal, AbsdK, yLeft, yRight, KLeft, KRight};

	F.function = &TrajectoryIntercept;
	F.params = &params;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc (T);
	gsl_root_fsolver_set (s, &F, x_lo, x_hi);

//	printf("\nUsing %s method:\n", gsl_root_fsolver_name (s));

//	printf ("%5s [%20s, %20s] %20s %21s %20s\n",
//		"iter", "lower", "upper", "root", 
//		"err", "err(est)");

	do{
		iter++;
		status = gsl_root_fsolver_iterate (s);
		r = gsl_root_fsolver_root (s);
		x_lo = gsl_root_fsolver_x_lower (s);
		x_hi = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);

		if (status == GSL_SUCCESS)
		printf ("Converged.\n");
//		printf ("Converged:\n");

//		printf ("%5d [%.18f, %.18f] %.18f %+.18f %.18f\n",
//			iter, x_lo, x_hi,
//			r, r - r_expected,
//			x_hi - x_lo);
	}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free (s);
//	printf ("\n");

	return r;
}



double
TrajectoryIntercept (double x, void *params) {
	double m,b;
	struct TrajectoryIntercept_params *p 
	= (struct TrajectoryIntercept_params *) params;

	double Ccval    = p->Ccval;
	double absdK    = p->absdK;
	double yleft    = p->yleft;
	double yright   = p->yright;
	double Kleft    = p->Kleft;
	double Kright   = p->Kright;

	m = (yright-yleft)/absdK;
	b = (yleft*Kright-yright*Kleft)/absdK;

	return exp(-PISQ*Ccval*(Kleft+x)) - m*(Kleft+x) - b;
}




//  Program Notes  //
//=========================================================================//
/*

See further documentation (UCLAresearchNotes.pdf, UCLAresearchPrograms.pdf, vlt_ThermStates.c, ring1.c) for elaboration.

*/
