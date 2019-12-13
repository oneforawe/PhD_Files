//  File Comments  //
//=========================================================================//

// AHHHHHHHHHHHHHHHHHHHHHH!!!!!!!!!!!!!!!!!!!!!
// I just took the derivative incorrectly.  Emin Menachekanian acted as my "research TA" and help me to catch my mistake.


/* FILENAME: EquilCompare_2DKT_FPquench.c */
/* VERSION: 2 (2011 Apr 06 - ...)
            Gamma should no longer fall below its new equilibrium curve */
/* AUTHOR: Andrew Forrester <aforrester@ucla.edu> */
/* DESCRIPTION:
   * This C program calculates time-dependent properties of a 2D system of superfluid (a thin layer of liquid 4He) before and after an instantaneous quench (temperature drop) starting at temperature Ti (Ti/Tkt = iTfrac) at time t=0 and going to temperature Tq (Tq/Tkt = qTfrac) at the next time step.  The system starts in equilibrium and the quench puts the system out of equilibrium so the properties evolve toward a new equilibrium.
   * The properties calculated are superfluid ratios (K, K/K0 = rho/rho0) and vortex-pair probability density G at various length-scales (l).
*/
/* EXT FILES: none */
/* COMPILE NOTES:
   * To compile, type "g++ -lm 2DKT_FPquench_adapt.c" without the quotes; then to run, type "./a.out".
   * NOTE: Be sure...
*/



//  Function Preparation  //
//=========================================================================//

// Standard routine header files //
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

// Constants //
const double PI   = 3.14159265358979323846;
const double PISQ = 9.86960440108935799230;
const double PICU = 31.00627668029981620634;
const double B    = 4.0*PICU;

// Program parameters/inputs //
const double lmax     = 10.0;         //  1  //  maximum of the length scale l (potentially approximate)
const int    lsteps   = 1250;         //  1250, 5000, 100000 //  lsteps = lmax/dl = 10/0.0001 = 100,000  (make sure this is divisible by ldatamax)
const int    lpts     = lsteps+1;     //  from l=0 to l=lmax, inclusive (incl. the boundary-condition-enforcing pnt, there are lpts+1 "l" pnts)
const int    ldatamax = 250;          //  50,100 //  max number of data points recorded per time step
const double iTfrac   = 0.95;         //  initial temperature fraction Ti/Tkt
const double a0       = 1.0;          //  a0 in units of a0 (!)
const double a04      = a0*a0*a0*a0;  //  a0 to the fourth power
const double K0c      = 0.747853;     //  critical (T=Tkt) bare superfluid ratio (or "coupling constant")

// Global variables //
int i;
double K[lpts+1], G[lpts+1];
double dGdt[1];



//  Function prototypes  //
//=========================================================================//

// rk4((*f)(),*x,y[],h,n);
// EqRecRel(x,z,dzdx,n);
// EquilMine(x,z,dzdx,n);
// EquilChus(x,z,dzdx,n);
void rk4(void (*f)(double, double*, double*, unsigned int), double *x, double y[], double h, unsigned n);
void EqRecRel(double x, double z[2], double dzdx[2], unsigned n);
void EquilMine(double x, double z[2], double dzdx[2], unsigned n);
void EquilChus(double x, double z[2], double dzdx[2], unsigned n);



//  Function definitions  //
//=========================================================================//

main(){
	// Main function definitions //
	int64_t n;
	int j,k;
	double dblsteps=lsteps, dl=lmax/dblsteps;
	double x, l[lpts];
	double K0, y0, G0, z[2];
	double Kcalc[lpts+1], Gcalc[lpts+1];
	double k1[lpts], k2[lpts], k3[lpts];
	double testG1[lpts+1], testG2[lpts+1];
	double t, dt, ttrig;
	double vpintegral, vpdens;
	FILE *outfile1;  // recording functions versus length-scale l ("vsl"), at selected times
	FILE *outfile2;  // recording functions versus time t ("vst"), at largest length scale lmax
	FILE *outfile1a;  // recording functions versus length-scale l ("vsl"), at selected times
	FILE *outfile2a;  // recording functions versus time t ("vst"), at largest length scale lmax
	FILE *outfile1b;  // recording functions versus length-scale l ("vsl"), at selected times
	FILE *outfile2b;  // recording functions versus time t ("vst"), at largest length scale lmax
	char *filename1, *filename2;
	char *filename1a, *filename2a;
	char *filename1b, *filename2b;

	// Prepare output files, print identification and values //
	// original:
	asprintf(&filename1, "EquilCompare_2DKT_FPquench_orig_vsl_T_%g_lmax%g_dl%g.dat", iTfrac, lmax, dl);
	asprintf(&filename2, "EquilCompare_2DKT_FPquench_orig_vst_T_%g_lmax%g_dl%g.dat", iTfrac, lmax, dl);
	outfile1 = fopen(filename1,"w");  //  E.g., "2DKT_FPquench.out"
	outfile2 = fopen(filename2,"w");  //  E.g., "2DKT_FPquench.out"
	fprintf(outfile1,"# Filename: %s\n", filename1);
	fprintf(outfile2,"# Filename: %s\n", filename2);
	fprintf(outfile1,"# Source: EquilCompare_2DKT_FPquench.c\n");
	fprintf(outfile2,"# Source: EquilCompare_2DKT_FPquench.c\n");
	fprintf(outfile1,"# Source version: %s\n", "2 (2011 Apr 06 - ...)");
	fprintf(outfile2,"# Source version: %s\n", "2 (2011 Apr 06 - ...)");
	fprintf(outfile1,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g, B=%21.21g\n", PI,PISQ,PICU,B);
	fprintf(outfile2,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g, B=%21.21g\n", PI,PISQ,PICU,B);
	fprintf(outfile1,"# Parameter values: iTfrac=%g, lmax=%g, lsteps=%i, (dl=%g), lpts=%i, ldatamax=%i, a0=%g, a04=%g, K0c=%g\n", iTfrac,lmax,lsteps,dl,lpts,ldatamax,a0,a04,K0c);
	fprintf(outfile2,"# Parameter values: iTfrac=%g, lmax=%g, lsteps=%i, (dl=%g), lpts=%i, ldatamax=%i, a0=%g, a04=%g, K0c=%g\n", iTfrac,lmax,lsteps,dl,lpts,ldatamax,a0,a04,K0c);

	// mine:
	asprintf(&filename1a, "EquilCompare_2DKT_FPquench_mine_vsl_T_%g_lmax%g_dl%g.dat", iTfrac, lmax, dl);
	asprintf(&filename2a, "EquilCompare_2DKT_FPquench_mine_vst_T_%g_lmax%g_dl%g.dat", iTfrac, lmax, dl);
	outfile1a = fopen(filename1a,"w");  //  E.g., "2DKT_FPquench.out"
	outfile2a = fopen(filename2a,"w");  //  E.g., "2DKT_FPquench.out"
	fprintf(outfile1,"# Filename: %s\n", filename1a);
	fprintf(outfile2,"# Filename: %s\n", filename2a);
	fprintf(outfile1,"# Source: EquilCompare_2DKT_FPquench.c\n");
	fprintf(outfile2,"# Source: EquilCompare_2DKT_FPquench.c\n");
	fprintf(outfile1,"# Source version: %s\n", "2 (2011 Apr 06 - ...)");
	fprintf(outfile2,"# Source version: %s\n", "2 (2011 Apr 06 - ...)");
	fprintf(outfile1,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g, B=%21.21g\n", PI,PISQ,PICU,B);
	fprintf(outfile2,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g, B=%21.21g\n", PI,PISQ,PICU,B);
	fprintf(outfile1,"# Parameter values: iTfrac=%g, lmax=%g, lsteps=%i, (dl=%g), lpts=%i, ldatamax=%i, a0=%g, a04=%g, K0c=%g\n", iTfrac,lmax,lsteps,dl,lpts,ldatamax,a0,a04,K0c);
	fprintf(outfile2,"# Parameter values: iTfrac=%g, lmax=%g, lsteps=%i, (dl=%g), lpts=%i, ldatamax=%i, a0=%g, a04=%g, K0c=%g\n", iTfrac,lmax,lsteps,dl,lpts,ldatamax,a0,a04,K0c);

	// Chu's:
	asprintf(&filename1b, "EquilCompare_2DKT_FPquench_chus_vsl_T_%g_lmax%g_dl%g.dat", iTfrac, lmax, dl);
	asprintf(&filename2b, "EquilCompare_2DKT_FPquench_chus_vst_T_%g_lmax%g_dl%g.dat", iTfrac, lmax, dl);
	outfile1b = fopen(filename1b,"w");  //  E.g., "2DKT_FPquench.out"
	outfile2b = fopen(filename2b,"w");  //  E.g., "2DKT_FPquench.out"
	fprintf(outfile1,"# Filename: %s\n", filename1b);
	fprintf(outfile2,"# Filename: %s\n", filename2b);
	fprintf(outfile1,"# Source: EquilCompare_2DKT_FPquench.c\n");
	fprintf(outfile2,"# Source: EquilCompare_2DKT_FPquench.c\n");
	fprintf(outfile1,"# Source version: %s\n", "2 (2011 Apr 06 - ...)");
	fprintf(outfile2,"# Source version: %s\n", "2 (2011 Apr 06 - ...)");
	fprintf(outfile1,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g, B=%21.21g\n", PI,PISQ,PICU,B);
	fprintf(outfile2,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g, B=%21.21g\n", PI,PISQ,PICU,B);
	fprintf(outfile1,"# Parameter values: iTfrac=%g, lmax=%g, lsteps=%i, (dl=%g), lpts=%i, ldatamax=%i, a0=%g, a04=%g, K0c=%g\n", iTfrac,lmax,lsteps,dl,lpts,ldatamax,a0,a04,K0c);
	fprintf(outfile2,"# Parameter values: iTfrac=%g, lmax=%g, lsteps=%i, (dl=%g), lpts=%i, ldatamax=%i, a0=%g, a04=%g, K0c=%g\n", iTfrac,lmax,lsteps,dl,lpts,ldatamax,a0,a04,K0c);

	// Boundary condition enforcement //
	K[lpts] = G[lpts] = 0.0;  // K=G=0 at l=lmax+dl, terminates recursion relations at l=lmax


	//////////////////////////////////////////////
	//                                          //
	//  (n=0) Initial Equilibrium Calculations  //
	//                                          //
	//////////////////////////////////////////////

	// Initialize time quantities //
	n = 0;
	t = 0.0;

	// Initialize quantities at smallest length-scale (l=0) //
	x = l[0] = 0.0;
	z[0] = K[0] = K0 = K0c/iTfrac;
	z[1] = y0 = exp(-PISQ*K0/2.0);
	G[0] = G0 = z[1]*z[1]*exp(-4.0*l[0])/a04;  // G (not = 0.0)
	vpintegral = 0.0;
	vpintegral += G[0]*exp(2.0*l[0])*dl;

	// Print initial data to screen and output files //
	printf(          "\n# time step n = %ld\tt = %e\t(dt = %e)\n", n,t,dt);
	fprintf(outfile1,"\n# time step n = %ld\tt = %e\t(dt = %e)\n", n,t,dt);
	printf(          "%s\t%s\t%s\t%s\t%s\t%s\n",   "T/Tc","1-T/Tc","l","K","K/K0","G");
	fprintf(outfile1,"# %s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","G");
	fprintf(outfile2,"# %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "time-step", "t", "T/Tc","1-T/Tc","l","K","K/K0","G","vpdens");
	printf(          "%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[0], K[0], K[0]/K0, G[0]);
	fprintf(outfile1,"%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[0], K[0], K[0]/K0, G[0]);

	// Calculate state-variables at progressively larger length scales, up to coherence length, via Runge-Kutta method (rk4) //
	for(i=1; i<lpts; i++){
		// Step out in length scale l, calculating next K, y, l, and G //
		rk4(EqRecRel, &x, z, dl, 2);  // equil: rk4 2D (K and y), EqRecRel
		l[i] = l[i-1] + dl;
		K[i] = z[0]; 
		G[i] = z[1]*z[1]*exp(-4.0*l[i])/a04;
		vpintegral += G[i]*exp(2.0*l[i])*dl;
		// Print at most [ldatamax] (e.g., 50) data points to screen and output file 1 //
		if(i%(lsteps/ldatamax)==0){
			printf(          "%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, G[i]);
			fprintf(outfile1,"%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, G[i]);
		}
	}
	vpdens = 2.0*PI*a0*a0*vpintegral;
	fflush(outfile1);

	// Print out the data at largest length scale to file 2 //
	fprintf(outfile2,"%ld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", n, t, iTfrac, 1-iTfrac, l[lpts-1], K[lpts-1], K[lpts-1]/K0, G[lpts-1], vpdens);
	fflush(outfile2);

	//////////////////////////////////////////////
	//                                          //
	//  Equilibrium Calculations (My equil eqn) //
	//                                          //
	//////////////////////////////////////////////

	printf("# My equilibrium equation\n");
	// Initialize quantities at smallest length-scale (l=0) //
	x = l[0] = 0.0;
	z[0] = K[0] = K0 = K0c/iTfrac;
	z[1] = G[0] = G0 = exp(-PISQ*K0)/a04;
	vpintegral = 0.0;
	vpintegral += G[0]*exp(2.0*l[0])*dl;

	// Print initial data to screen and output files //
	printf(           "%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[0], K[0], K[0]/K0, G[0]);
	fprintf(outfile1a,"%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[0], K[0], K[0]/K0, G[0]);

	// Calculate state-variables at progressively larger length scales, up to coherence length, via Runge-Kutta method (rk4) //
	for(i=1; i<lpts; i++){
		// Step out in length scale l, calculating next K, y, l, and G //
		rk4(EquilMine, &x, z, dl, 2);  // equil: rk4 2D (K and y), EqRecRel
		l[i] = l[i-1] + dl;
		K[i] = z[0]; 
		G[i] = z[1];
		vpintegral += G[i]*exp(2.0*l[i])*dl;
		// Print at most [ldatamax] (e.g., 50) data points to screen and output file 1 //
		if(i%(lsteps/ldatamax)==0){
			printf(           "%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, G[i]);
			fprintf(outfile1a,"%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, G[i]);
		}
	}
	vpdens = 2.0*PI*a0*a0*vpintegral;
	fflush(outfile1a);

	// Print out the data at largest length scale to file 2 //
	fprintf(outfile2a,"%ld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", n, t, iTfrac, 1-iTfrac, l[lpts-1], K[lpts-1], K[lpts-1]/K0, G[lpts-1], vpdens);
	fflush(outfile2a);

	//////////////////////////////////////////////
	//                                          //
	//  (n=0) Initial Equilibrium (Chu eq eqn)  //
	//                                          //
	//////////////////////////////////////////////

	printf("# Chu's equilibrium equation\n");
	// Initialize quantities at smallest length-scale (l=0) //
	x = l[0] = 0.0;
	z[0] = K[0] = K0 = K0c/iTfrac;
	z[1] = G[0] = G0 = exp(-PISQ*K0)/a04;
	vpintegral = 0.0;
	vpintegral += G[0]*exp(2.0*l[0])*dl;

	// Print initial data to screen and output files //
	printf(           "%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[0], K[0], K[0]/K0, G[0]);
	fprintf(outfile1b,"%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[0], K[0], K[0]/K0, G[0]);

	// Calculate state-variables at progressively larger length scales, up to coherence length, via Runge-Kutta method (rk4) //
	for(i=1; i<lpts; i++){
		// Step out in length scale l, calculating next K, y, l, and G //
		rk4(EquilChus, &x, z, dl, 2);  // equil: rk4 2D (K and y), EqRecRel
		l[i] = l[i-1] + dl;
		K[i] = z[0]; 
		G[i] = z[1];
		vpintegral += G[i]*exp(2.0*l[i])*dl;
		// Print at most [ldatamax] (e.g., 50) data points to screen and output file 1 //
		if(i%(lsteps/ldatamax)==0){
			printf(           "%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, G[i]);
			fprintf(outfile1b,"%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, G[i]);
		}
	}
	vpdens = 2.0*PI*a0*a0*vpintegral;
	fflush(outfile1b);

	// Print out the data at largest length scale to file 2 //
	fprintf(outfile2b,"%ld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", n, t, iTfrac, 1-iTfrac, l[lpts-1], K[lpts-1], K[lpts-1]/K0, G[lpts-1], vpdens);
	fflush(outfile2b);



	fclose(outfile1);
	fclose(outfile2);
	fclose(outfile1a);
	fclose(outfile2a);
	fclose(outfile1b);
	fclose(outfile2b);
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
	for(j=0;j<n;j++)  y[j] += h*(k[0][j]+2*(k[1][j]+k[2][j])+k[3][j])/6;
	*x += h;
	for(i=0;i<4;i++)  free((char *) k[i]);
	free((char *)s);
}



// equilibrium (K,y) recursion relations //
void EqRecRel(double x, double z[2], double dzdx[2], unsigned n){
	dzdx[0] = -B*z[1]*z[1]*z[0]*z[0];
	dzdx[1] = (2-PI*z[0])*z[1];
}



// My equilibrium (K,G) recursion relations //
void EquilMine(double x, double z[2], double dzdx[2], unsigned n){
	dzdx[0] = -a04*B*z[1]*z[0]*z[0]*exp(4.0*x);
	dzdx[1] = (4.0-4.0*x-2.0*PI*z[0])*z[1];
}



// Chu's equilibrium (K,G) recursion relations //
void EquilChus(double x, double z[2], double dzdx[2], unsigned n){
	dzdx[0] = -a04*B*z[1]*z[0]*z[0]*exp(4.0*x);
	dzdx[1] = -2.0*PI*z[0]*z[1];
}



//  Program Notes  //
//=========================================================================//
/*

See further documentation (UCLAresearchNotes.pdf, UCLAresearchPrograms.pdf, derk1.c) for elaboration.



vortex pair density
vpdens = int_a0^amax G*d^2a = int_a0^amax int_0^2PI G*a*dtheta*da = 2PI int_a0^amax G*a*da = 
a = a0*exp(l)
da = a0*exp(l)*dl
a*da = a0^2*exp(2l)*dl

vpdens = 2PI*a0^2 int_0^lmax G*exp(2l)*dl
vpintegral += G*exp(2l)*dl
vpdens = 2PI*a0^2*vpintegral

*/
