//  File Comments  //
//=========================================================================//

/* FILENAME: 2DKT_FPquench_adapt.c */
/* VERSION: 1 (2011 Jan 02)
/* AUTHOR: Andrew Forrester <aforrester@ucla.edu> */
/* DESCRIPTION:
   * This C program calculates time-dependent properties of a 2D system of superfluid (a thin layer of liquid 4He) before and after an instantaneous quench (temperature drop) starting at temperature Ti (Ti/Tkt = iTfrac) at time t=0 and going to temperature Tq (Tq/Tkt = qTfrac) at the next time step.  The system starts in equilibrium and the quench puts the system out of equilibrium so the properties evolve toward a new equilibrium.
   * The properties calculated are superfluid ratios (K, K/K0 = rho/rho0) and vortex-pair probability density G at various length-scales (l).
*/
/* EXT FILES: none */
/* COMPILE NOTES:
   * To compile, type "g++ -lm 2DKT_FPquench_adapt.c" without the quotes; then to run, type "./a.out".
   * If you want to have more info for debugging, using dgb,
     to compile, type "g++ -lm -g 2DKT_FPquench_adapt.c" without the quotes;
     then to run and get information about how long it took to run, type "time ./a.out".
   
   * To compile, type (without the quotes) "g++ -Wall -I/usr/local/src/libslack-0.6 -c 2DKT_FPquench_adapt.c"       works?
   * To link to the gsl libraries, type    "g++ -g -L/usr/lib -lm -lslack 2DKT_FPquench_adapt.o"                    works?
   * To link to the gsl libraries, type    "g++ -g -L/usr/local/src/libslack-0.6 -lm -lslack 2DKT_FPquench_adapt.o" works?
   * To link to the gsl libraries, type    "g++ -g -L/usr/lib -lm -lgsl -lgslcblas 2DKT_FPquench_adapt.o"
   * To run, type                          "time ./a.out".

   * NOTE: In GNU/Linux, you may have to make the stack size unlimited (to allow lsteps to be a large integer, approx >150,000):
     Check limits with the command "ulimit -a", and set the stack size to unlimited with "ulimit -s unlimited"
*/
/* POSSIBLE IMPROVEMENTS:
   * The error could be accumulated through each time step, not just through each l-step, which is already done.  The error at each l-point would have to be recorded and added with each time-step.  To be nit-picky, we should also have a formula for the error at l=0.  (For now I set the error at l=0 to be zero.)
*/



//  Function Preparation  //
//=========================================================================//

// Standard routine header files //
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
//#include <sys/resource.h>

// Constants //
const double PI   = 3.14159265358979323846;
const double PISQ = 9.86960440108935799230;
const double PICU = 31.00627668029981620634;
const double B    = 4.0*PICU;

// Program parameters/inputs //
const double lmax      = 10.0;      //  1  //  maximum of the length scale l (potentially approximate)
const int    lsteps    = 500000;    //  5000 // 100000 //  lsteps = lmax/dl0
                                    //  the higher the accuracy you want (the lower DRE and TOL) the more steps you'll have to take => larger lsteps
                                    //  lmax/dl0 = 10.0/0.0001 = 100,000
                                    //  lmax/dl0 = 10.0/0.008  =   1,250
                                    //  lmax/dl0 = 10.0/(2e-5) = 500,000
const int    lpts      = lsteps+1;  //  from l=0 to l=lmax, inclusive (incl. the boundary-condition-enforcing pnt, there are lpts+1 "l" pnts)
const int    ldatamax  = 50;        //  maximum number of l data points recorded per time step (beyond the l=0 data point, not including error bar "min" and "max" points)
const double iTfrac    = 0.95;      //  initial temperature fraction Ti/Tkt
const double qTfrac    = 0.10;      //  quench temperature fraction Tq/Tkt
const double tmax      = 10000.0;   //  20600 //  max unitless time (see if(t==...) below)
const double dt0       = 1.8e-4;    //  the time increment (in units of the "diffusion time")
const double a0        = 1.0;       //  a0 in units of a0 (!)
const double a04       = 1.0;       //  a0 to the fourth power
const double K0c       = 0.747853;  //  critical (T=Tkt) bare superfluid ratio (or "coupling constant")
const double DRE       = 1e-15;     //  desired relative error
                         //  I GET A SEGMENTATION FAULT IF I USE  DRE = 1e-16
                         //  BUT THE TOTAL ERROR (IN G) IS BIGGER THAN THE (G) VALUE
const double TOL       = 1e-5;//1e-26;
                         //  IT SEEMS TO TAKE FOREVER WITH TOL = 1e-27
const double PSHRINK   = -0.25;
const double PGROW     = -0.20;
const double MAXFACTOR = 3.0;
const double SAFETY    = 0.9;
const int    FILENAME_BUFFER_SIZE = 200;

// rkf45: adaptive-step-size Runge-Kutta-Fehlberg method, 5th order
const double b21=0.25,
	     b31=3.0/32.0,      b32=9.0/32.0,
	     b41=1932.0/2197.0, b42=-7200.0/2197.0, b43=7296.0/2197.0,
	     b51=439.0/216.0,   b52=-8.0,           b53=3680.0/513.0,     b54=-845.0/4104.0,
	     b61=-8.0/27.0,     b62=2.0,            b63=-3544.0/2565.0,   b64=1859.0/4104.0,    b65=-11.0/40.0,
	      a1=0.0,            a2=0.25,            a3=0.375,             a4=12.0/13.0,         a5=1.0,         a6=0.5,
	      c1=16.0/135.0,     c2=0.0,             c3=6656.0/12825.0,    c4=28561.0/56430.0,   c5=-9.0/50.0,   c6=2.0/55.0,
	     dc1=c1-25.0/216.0, dc2=c2-0.0,         dc3=c3-1408.0/2565.0, dc4=c4-2197.0/4104.0, dc5=c5+1.0/5.0, dc6=c6-0.0;

/*
// rkck45: adaptive-step-size Runge-Kutta-Cash-Karp method, 5th order
const double b21=0.2,
	     b31=3.0/40.0,          b32=9.0/40.0,
	     b41=0.3,               b42=-0.9,        b43=1.2,
	     b51=-11.0/54.0,        b52=2.5,         b53=-70.0/27.0,         b54=35.0/27.0,
	     b61=1631.0/55296.0,    b62=175.0/512.0, b63=575.0/13824.0,      b64=44275.0/110592.0,   b65=253.0/4096.0,
	      a1=0.0,                a2=0.2,          a3=0.3,                 a4=0.6,                 a5=1.0,              a6=0.875,
	      c1=37.0/378.0,         c2=0.0,          c3=250.0/621.0,         c4=125.0/594.0,         c5=0.0,              c6=512.0/1771.0,
	     dc1=c1-2825.0/27648.0, dc2=0.0,         dc3=c3-18575.0/48384.0, dc4=c4-13525.0/55296.0, dc5=-277.00/14336.0, dc6=c6-0.25;
*/

//struct rlimit{
//	rlim_t rlim_cur;  // Soft limit
//	rlim_t rlim_max;  // Hard limit (ceiling for rlim_cur)
//};

// Global variables //
int i;
double l[lpts], dl[lpts], K[lpts+1], G[lpts+1];
double dGdt[1];



//  Function prototypes  //
//=========================================================================//

// Ky_rkf45((*derivs)(),*x,h,y[],n,dre,*htemp,ytemp[],yerr[],*maxrelerr);
// K_rkf45_NoAdjust((*derivs)(),*x,h,y[],n,dre,yerr[],*maxrelerr);
// G_rkf45(imax,l[lpts],dl[lpts],G[lpts+1],K[lpts+1],*t,*dt,Gerr[lpts]);
// EqRecRel(x,z[2],dzdx[2],n);
// fpKRecRel(x,Kcalc[1],dKdx[1],n);
// fpGRecRel(x,dx1,dx2,Gcalc[lpts+1],K[lpts+1],*dGdt);
void Ky_rkf45( void (*derivs)(double, double*, double*, unsigned),
	       double *x, double h, double y[], unsigned n, double dre, double *htemp, double ytemp[], double yerr[], double *maxrelerr );
void K_rkf45_NoAdjust( void (*derivs)(double, double*, double*, unsigned),
		       double *x, double h, double y[], unsigned n, double dre, double yerr[], double *maxrelerr );
void G_rkf45(int imax, double l[lpts], double dl[lpts], double G[lpts+1], double K[lpts+1], double *t, double *dt, double Gerr[lpts]);
void EqRecRel(double x, double z[2], double dzdx[2], unsigned n);
void fpKRecRel(double x, double Kcalc[1], double dKdx[1], unsigned n);
void fpGRecRel(double x, double dx1,  double dx2, double Gcalc[lpts+1], double K[lpts+1], double *dGdt);

//int getrlimit(int resource, struct rlimit *rlim);
//int setrlimit(int resource, const struct rlimit *rlim);



//  Function definitions  //
//=========================================================================//

int main(){
	// In GNU/Linux, you may have to do this to allow lsteps to be a large integer (approx >150,000)  //
	//struct rlimit rlimit_stack;
	//getrlimit(RLIMIT_STACK,&rlimit_stack);
	//rlimit_stack.rlim_max = RLIM_INFINITY;
	//rlimit_stack.rlim_cur = RLIM_INFINITY;
	//setrlimit(RLIMIT_STACK,&rlimit_stack);
	// You can find this information with the command "man setrlimit"
	// Or, an alternative to this block of code is to check limits with the command "ulimit -a", and set the stack size to unlimited with "ulimit -s unlimited"

	// Main function definitions //
	int64_t n;
	int imax,accuracy,j,k;
	double t, dt;
	double dblsteps=lsteps, dbldatamax=ldatamax;
	double dl0=lmax/dblsteps, dataDl=lmax/dbldatamax;
	double dltemp;
	double K0, y0;
	double z[2], ztemp[2], zerr[2];
	double Kcalc[1], Kcalctemp[1], Kcalcerr[1];
	double Gcalc[1], Gcalctemp[1], Gcalcerr[1];
	double Ktemp[lpts], Kerr[lpts];
	double Gtemp[lpts], Gerr[lpts];
	double KerrTot, GerrTot;
	double maxrelerror, Grelerror, factor_dl;
//	double k1[lpts], k2[lpts], k3[lpts];
//	double testG1[lpts+1], testG2[lpts+1];
	FILE *outfile;
	//char filename[FILENAME_BUFFER_SIZE];
	//int filename_num;
	char *filename;

	// Prepare output file, print identification and values //
	if( asprintf(&filename, "2DKT_FPquench_adapt_lmax%g_dl0%g_tmax%g_dt0%g_DRE%g_TOL%g.dat", lmax, dl0, tmax, dt0, DRE, TOL) < 0 ){
		printf("%s","Filename error.");
		return EXIT_FAILURE;
	}
	outfile = fopen(filename,"w");  //  E.g., "2DKT_FPquench_adapt.out"
	fprintf(outfile,"# Filename: %s\n", filename);
	fprintf(outfile,"# Source: 2DKT_FPquench_adapt.c\n");
	fprintf(outfile,"# Source version: %s\n", "1 (2011 Jan 02)");
	fprintf(outfile,"# Parameter values: lmax=%g, lsteps=%i, (dl0=%g), lpts=%i, ldatamax=%i, iTfrac=%g, qTfrac=%g, tmax=%g, dt0=%g, a0=%g, a04=%g, K0c=%g, DRE=%g, TOL=%g, PSHRINK=%g, PGROW=%g, MAXFACTOR=%g, SAFETY=%g\n", lmax,lsteps,dl0,lpts,ldatamax,iTfrac,qTfrac,tmax,dt0,a0,a04,K0c,DRE,TOL,PSHRINK,PGROW,MAXFACTOR,SAFETY);

	// Prepare output file, print identification and values //
	//filename_num = snprintf(filename, FILENAME_BUFFER_SIZE, "2DKT_FPquench_adapt_lmax%g_dl0%g_tmax%g_dt0%g_DRE%g_TOL%g.dat", lmax, dl0, tmax, dt0, DRE, TOL);
	//if(filename_num>=FILENAME_BUFFER_SIZE){
	//	printf("%s","Filename size error.");
	//	return EXIT_FAILURE;
	//}
	//outfile = fopen(filename,"w");  //  E.g., "2DKT_FPquench_adapt.out"
	//fprintf(outfile,"# Filename: %s\n", filename);
	//fprintf(outfile,"# Source: 2DKT_FPquench_adapt.c\n");
	//fprintf(outfile,"# Source version: %s\n", "1 (2011 Jan 02)");
	//fprintf(outfile,"# Parameter values: lmax=%g, lsteps=%i, (dl0=%g), lpts=%i, ldatamax=%i, iTfrac=%g, qTfrac=%g, tmax=%g, dt0=%g, a0=%g, a04=%g, K0c=%g, DRE=%g, TOL=%g, PSHRINK=%g, PGROW=%g, MAXFACTOR=%g, SAFETY=%g\n", lmax,lsteps,dl0,lpts,ldatamax,iTfrac,qTfrac,tmax,dt0,a0,a04,K0c,DRE,TOL,PSHRINK,PGROW,MAXFACTOR,SAFETY);

	// Prepare output file, print identification and values //
	//sprintf(filename, "2DKT_FPquench_adapt_lmax%g_dl0%g_tmax%g_dt0%g_DRE%g_TOL%g.dat", lmax, dl0, tmax, dt0, DRE, TOL);
	//outfile = fopen(filename,"w");  //  E.g., "2DKT_FPquench_adapt.out"
	//fprintf(outfile,"# Filename: %s\n", filename);
	//fprintf(outfile,"# Source: 2DKT_FPquench_adapt.c\n");
	//fprintf(outfile,"# Source version: %s\n", "1 (2011 Jan 02)");
	//fprintf(outfile,"# Parameter values: lmax=%g, lsteps=%i, (dl0=%g), lpts=%i, ldatamax=%i, iTfrac=%g, qTfrac=%g, tmax=%g, dt0=%g, a0=%g, a04=%g, K0c=%g, DRE=%g, TOL=%g, PSHRINK=%g, PGROW=%g, MAXFACTOR=%g, SAFETY=%g\n", lmax,lsteps,dl0,lpts,ldatamax,iTfrac,qTfrac,tmax,dt0,a0,a04,K0c,DRE,TOL,PSHRINK,PGROW,MAXFACTOR,SAFETY);


	//////////////////////////////////////////////
	//                                          //
	//  (n=0) Initial Equilibrium Calculations  //
	//                                          //
	//////////////////////////////////////////////

	// Initialize time quantities //
	n = 0;
	t = 0.0;
	dt = dt0;  // We won't step in time until the quench.

	// Initialize quantities at smallest length-scale (l=0) //
	l[0] = 0.0;
	dl[0] = dl0;
	K0 = K0c/iTfrac;
	y0 = exp(-PISQ*K0/2);
	z[0] = K0;
	z[1] = y0;
	K[0] = z[0];
	G[0] = z[1]*z[1]*exp(-4.0*l[0])/a04;  // G (not = 0.0)
	Kerr[0] = 0;  // There should be some kind of thermodynamic formula for this error.
	Gerr[0] = 0;  // There should be some kind of thermodynamic formula for this error.

	// Print headings and initial data to screen and output file //
	printf(         "\n# time step n = %ld\tt = %e\t(dt = %e)\n", n,t,dt);
	fprintf(outfile,"\n# time step n = %ld\tt = %e\t(dt = %e)\n", n,t,dt);
	printf(         "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",  "Data-label", "lstep", "T/Tc","1-T/Tc","l","K","K/K0","G");
	fprintf(outfile,"# %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",  "Data-label", "lstep", "T/Tc","1-T/Tc","l","K","K/K0","G");
	printf(         "%i\t%i\t%g\t%g\t%g\t%g\t%g\t%g\n", 0, 0, iTfrac, 1-iTfrac, l[0], K[0], K[0]/K0, G[0]);
	fprintf(outfile,"%i\t%i\t%g\t%g\t%g\t%g\t%g\t%g\n", 0, 0, iTfrac, 1-iTfrac, l[0], K[0], K[0]/K0, G[0]);

	// Calculate state-variables at progressively larger length scales, up to coherence length, via Runge-Kutta method //
	i = 0; // for stepping l (i = "lstep")
	j = 1; // for selecting data to print to file
	k = 1; // for data labelling
	while(l[i]<lmax){ // && i<2){
		// Increment to next lstep, with unknown values for l, K, y, G //
		i++;
////		printf("\ni = %i\n\n",i);
		// Step out in length scale l, calculating next K, y, l, G (and errors in K and G) //
		accuracy = 0;
		while(accuracy==0){
			Ky_rkf45(EqRecRel, &l[i-1], dl[i-1], z, 2, DRE, &dltemp, ztemp, zerr, &maxrelerror);
			// equil: rkf6 2D (K and y), EqRecRel, adaptive step to desired relative error DRE
			//printf("%s\n","made it out of rk!");
			//printf("\tl[%i-1] = %g\n",i,l[i-1]);
			//printf("\tdltemp = %g\n",dltemp);
			l[i] = l[i-1]+dltemp;
			//printf("\tl[%i] = %g\n",i,l[i]);
			//printf("\tztemp[1] = %g\n",ztemp[1]);
			//printf("\texp(-4.0*l[%i]) = %g\n",i,exp(-4.0*l[i]));
			//printf("\tztemp[1]*ztemp[1]*exp(-4.0*l[%i])/a04 = %g\n",i,ztemp[1]*ztemp[1]*exp(-4.0*l[i])/a04);
			G[i] = ztemp[1]*ztemp[1]*exp(-4.0*l[i])/a04;
			//printf("\tG[%i] = %g\n",i,G[i]);
			Gerr[i] = zerr[1]*2*ztemp[1]*exp(-4.0*l[i])/a04;
			//printf("\tGerr[%i] = %g\n",i,Gerr[i]);
			// check error (in G), redo (w/ smaller increment) or proceed (w/o changing increment) //
			Grelerror = fabs( Gerr[i]/TOL );  // (calculated error)/(desired error)   ( Earlier, TOL->(DRE*G[i]) )
			if(Grelerror>maxrelerror)  maxrelerror = Grelerror;
			//printf("\terrfrac = %g\n",errfrac);
			if(maxrelerror>1.0){  // calculated error is too big
////				printf("%s","accuracy = 0\n");
				accuracy = 0;
				dl[i-1] = SAFETY*dltemp*pow(maxrelerror,PSHRINK);
				// Failure print out //
				//printf("G accuracy Fail: %g\t%g\t%g\t%g\n", l[i-1], dl[i-1], ztemp[0], G[i]);
			}
			else{  // calculated error small enough (or overly small), so accept current dl and grow it for next step
				//printf("%s","accuracy = 1\n");
				accuracy = 1;
				dl[i-1] = dltemp;
				factor_dl = SAFETY*pow(maxrelerror,PGROW);
				if(factor_dl<MAXFACTOR)  dl[i] = dltemp*factor_dl;
				else                     dl[i] = dltemp*MAXFACTOR;
				if(l[i]+dl[i]>lmax)      dl[i] = lmax-l[i];
				z[0] = ztemp[0];
				K[i] = z[0];
				Kerr[i] = zerr[0];
				z[1] = ztemp[1];
				// Success print out //
////				printf("G accuracy Success: %g\t%g\t%g\t%g\n", l[i-1], dl[i-1], ztemp[0], G[i]);
			}
		}
		// Print at most [ldatamax] (e.g., 50) data points to screen and output file //
		if(l[i]>=j*dataDl){
			while(l[i]>=j*dataDl)  j++;
			printf(         "%i\t%i\t%g\t%g\t%g\t%g\t%g\t%g\n", k, i, iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, G[i]);
			fprintf(outfile,"%i\t%i\t%g\t%g\t%g\t%g\t%g\t%g\n", k, i, iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, G[i]);
			k++;
		}
	}
	imax = i;

	// Calculate and print total error and error bar min and max for the final values //
	for(i=0;i<=imax;i++){
		GerrTot += Gerr[i];
		KerrTot += Kerr[i];
	}
	fprintf(outfile,"%s\t%i\t%g\t%g\t%g\t%g\t%g\t%g\n", "min", imax, iTfrac, 1-iTfrac, l[imax], (K[imax]-KerrTot), (K[imax]-KerrTot)/K0, (G[imax]-GerrTot) );
	fprintf(outfile,"%s\t%i\t%g\t%g\t%g\t%g\t%g\t%g\n", "max", imax, iTfrac, 1-iTfrac, l[imax], (K[imax]+KerrTot), (K[imax]+KerrTot)/K0, (G[imax]+GerrTot) );
	fprintf(outfile,"# Total error in K = %g\n", KerrTot);
	fprintf(outfile,"# Total error in G = %g\n", GerrTot);
	fflush(outfile);

	// Boundary condition enforcement for next stage of program //
	K[imax+1] = G[imax+1] = 0.0;  // K=G=0 at l=lmax+dl, terminates recursion relations at l=lmax


	//////////////////////////////////////////////////////////////
	//                                                          //
	//  (n>=1) After Instantaneous "Quench" (Temperature Drop)  //
	//                                                          //
	//////////////////////////////////////////////////////////////

	//  UNLESS I USE SOME KIND OF INTERPOLATION, I MUST USE THE SAME dl[] IN EVERY SUBSEQUENT l-STEP SINCE dGdt DEPENDS ON THE OLD G

	// Initialize the changed quantities at smallest length-scale (l=0) //
	K0 = K0c/qTfrac;
	//y0 = exp(-PISQ*K0/2);
	K[0] = K0;
	//G[0] = y0*y0*exp(-4.0*l[0])/a04;  // G (not = 0.0)

	// Instantaneous change in K //
	Kcalc[0] = K0;
	for(i=1; i<=imax; i++){
		//printf("%s\n","hello1!");
		K_rkf45_NoAdjust(fpKRecRel, &l[i-1], dl[i-1], Kcalc, 1, DRE, Kcalcerr, &maxrelerror);
		K[i] = Kcalc[0];
		Kerr[i] = Kcalcerr[0];
		if(maxrelerror>1.0)  printf("Large error in K!\tmaxrelerror = %g\tKerr[%i] = %g\n",maxrelerror,i,Kerr[i]);
	}

	// Advance in time //
	while(t<tmax && n<6){
		// //
		n++;
		// Put the following bit of code wherever it makes the most sense //
		printf(         "\n# time step n = %ld\tt = %e\t(dt = %e)\n", n,t,dt);
		fprintf(outfile,"# time step n = %ld\tt = %e\t(dt = %e)\n", n,t,dt);

		// Calculate dGdt for all length scales, and update G //

		// advance t calculating next G (and Gerr)
		G_rkf45(imax, l, dl, G, K, &t, &dt, Gerr);  // do something with Gerr, if you want!

		// get next K (which depends on the new G)  and print the results
		Kcalc[0] = K0;
		for(i=1; i<=imax; i++){
			K_rkf45_NoAdjust(fpKRecRel, &l[i-1], dl[i-1], Kcalc, 1, DRE, Kcalcerr, &maxrelerror);
			K[i] = Kcalc[0];
			Kerr[i] = Kcalcerr[0];
			if(maxrelerror>1.0)  printf("\n%s\n","Bogus!  Large error in K!");

			// Print at most [ldatamax] (e.g., 50) data points to screen and output file //
			if(l[i]>=j*dataDl){
				while(l[i]>=j*dataDl)  j++;
				printf(         "%i\t%i\t%g\t%g\t%g\t%g\t%g\t%g\n", k, i, iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, G[i]);
				fprintf(outfile,"%i\t%i\t%g\t%g\t%g\t%g\t%g\t%g\n", k, i, iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, G[i]);
				k++;
			}
		}

		// Calculate and print total error, min, and max final values //
		for(i=0; i<=imax; i++){
			GerrTot += Gerr[i];
			KerrTot += Kerr[i];
		}
		fprintf(outfile,"%s\t%i\t%g\t%g\t%g\t%g\t%g\t%g\n", "min", imax, iTfrac, 1-iTfrac, l[imax], (K[imax]-KerrTot), (K[imax]-KerrTot)/K0, (G[imax]-GerrTot) );
		fprintf(outfile,"%s\t%i\t%g\t%g\t%g\t%g\t%g\t%g\n", "max", imax, iTfrac, 1-iTfrac, l[imax], (K[imax]+KerrTot), (K[imax]+KerrTot)/K0, (G[imax]+GerrTot) );
		fprintf(outfile,"# Total error in K = %g\n", KerrTot);
		fprintf(outfile,"# Total error in G = %g\n", GerrTot);
		fflush(outfile);

	}

/*
		// At certain times, print out the results //
		if( (t>=0.1&&k==1) || (t>=1&&k==2) || (t>=10&&k==3) || (t>=100&&k==4) || (t>=1000&&k==5) || (t>=10000&&k==6) ){
			fprintf(outfile,"\n# time step n = %ld\tt = %e\t(dt = %e)\n", n,t,dt);
			fprintf(outfile,"# %s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","G");
			for(i=0; i<lpts; i++){
				// Print at most [ldatamax] (e.g., 1000) data points to screen and output file //
				if(i%(lsteps/ldatamax)==0){
					fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\n", qTfrac, 1-qTfrac, l[i], K[i], K[i]/K0, G[i]);
					printf(         "%g\t%g\t%g\t%g\t%g\t%g\n", qTfrac, 1-qTfrac, l[i], K[i], K[i]/K0, G[i]);
				}
			}
			fflush(outfile);
			k++;
		}

		printf(         "%g\t%g\t%g\n", l[lpts-1], K[lpts-1]/K0, G[lpts-1]);
		//fprintf(outfile,"%g\t%g\t%g\n", l[lpts-1], KG[lpts-1].c[0]/K0, KG[lpts-1].c[1]);
		fflush(outfile);

	}

*/
	fclose(outfile);
	return EXIT_SUCCESS;
}



void Ky_rkf45(void (*derivs)(double, double*, double*, unsigned), double *x, double h, double y[], unsigned n, double dre, double *htemp, double ytemp[], double yerr[], double *maxrelerr){
	int i,accurate=0;
	double *k[6],*ytest,temprelerr;

	for(i=0;i<6;i++)  k[i] = (double *)calloc(n,sizeof(double));  // the k's are slopes, dy/dx
	ytest = (double *)calloc(n,sizeof(double));

	*htemp = h;
	while(accurate==0){
		*maxrelerr=1e-300;
		derivs(*x,y,k[0],n);                // get k1=k[0]
		for(i=0;i<n;i++)  ytest[i] = y[i]+*htemp*(b21*k[0][i]);
		derivs(*x+*htemp*a2,ytest,k[1],n);  // get k2=k[1]
		for(i=0;i<n;i++)  ytest[i] = y[i]+*htemp*(b31*k[0][i]+b32*k[1][i]);
		derivs(*x+*htemp*a3,ytest,k[2],n);  // get k3=k[2]
		for(i=0;i<n;i++)  ytest[i] = y[i]+*htemp*(b41*k[0][i]+b42*k[1][i]+b43*k[2][i]);
		derivs(*x+*htemp*a4,ytest,k[3],n);  // get k4=k[3]
		for(i=0;i<n;i++)  ytest[i] = y[i]+*htemp*(b51*k[0][i]+b52*k[1][i]+b53*k[2][i]+b54*k[3][i]);
		derivs(*x+*htemp*a5,ytest,k[4],n);  // get k5=k[4]
		for(i=0;i<n;i++)  ytest[i] = y[i]+*htemp*(b61*k[0][i]+b62*k[1][i]+b63*k[2][i]+b64*k[3][i]+b65*k[4][i]);
		derivs(*x+*htemp*a6,ytest,k[5],n);  // get k6=k[5]

		for(i=0;i<n;i++){
			ytemp[i] = y[i] + *htemp*(c1*k[0][i]+c2*k[1][i]+c3*k[2][i]+c4*k[3][i]+c5*k[4][i]+c6*k[5][i]);
			//printf("\tytemp[%i] = %g\t",i,ytemp[i]);
			yerr[i] = *htemp*(dc1*k[0][i]+dc2*k[1][i]+dc3*k[2][i]+dc4*k[3][i]+dc5*k[4][i]+dc6*k[5][i]);
			//printf("yerr[%i] = %g\t",i,yerr[i]);
			temprelerr = fabs( yerr[i]/(dre*ytemp[i]) );  // (calculated error)/(desired error)   ( Earlier, (dre*ytemp[i])->TOL )
			//printf("temperr = %g\n",temperr);
			if(temprelerr>*maxrelerr)  *maxrelerr = temprelerr;
		}
		//printf("\tmaxerr = %g\n",maxerr);
		if(*maxrelerr>1.0){  // calculated error is too big
			accurate = 0;
			// Failure print out //
			//printf("rkfail: %g\t%g\t%g\t%g\t%g\t%g\t%g\t%i\n", *x, *htemp, ytemp[0], yerr[0], ytemp[1], yerr[1], maxerr, accurate);
			*htemp *= SAFETY*pow(*maxrelerr,PSHRINK);
		}
		else{  // calculated error small enough (potentially there are further checks outside the program, before growing h)
			accurate = 1;
			// (THIS IS DONE OUTSIDE OF THE PROGRAM!) factor = SAFETY*pow(*maxrelerr,PGROW);
			// if(factor<MAXFACTOR)  *htemp *= factor;
			// else                  *htemp *= MAXFACTOR;
		}
	}

	for(i=0;i<6;i++)  free((char *)k[i]);
	free((char *)ytest);
}



void K_rkf45_NoAdjust(void (*derivs)(double, double*, double*, unsigned), double *x, double h, double y[], unsigned n, double dre, double yerr[], double *maxrelerr){
	int i;
	double *k[6],*ytest;

	for(i=0;i<6;i++)  k[i] = (double *)calloc(n,sizeof(double));  // the k's are slopes, dy/dx
	ytest = (double *)calloc(n,sizeof(double));

	derivs(*x,y,k[0],n);           // get k1=k[0]
	for(i=0;i<n;i++)  ytest[i] = y[i]+h*(b21*k[0][i]);
	derivs(*x+h*a2,ytest,k[1],n);  // get k2=k[1]
	for(i=0;i<n;i++)  ytest[i] = y[i]+h*(b31*k[0][i]+b32*k[1][i]);
	derivs(*x+h*a3,ytest,k[2],n);  // get k3=k[2]
	for(i=0;i<n;i++)  ytest[i] = y[i]+h*(b41*k[0][i]+b42*k[1][i]+b43*k[2][i]);
	derivs(*x+h*a4,ytest,k[3],n);  // get k4=k[3]
	for(i=0;i<n;i++)  ytest[i] = y[i]+h*(b51*k[0][i]+b52*k[1][i]+b53*k[2][i]+b54*k[3][i]);
	derivs(*x+h*a5,ytest,k[4],n);  // get k5=k[4]
	for(i=0;i<n;i++)  ytest[i] = y[i]+h*(b61*k[0][i]+b62*k[1][i]+b63*k[2][i]+b64*k[3][i]+b65*k[4][i]);
	derivs(*x+h*a6,ytest,k[5],n);  // get k6=k[5]

	for(i=0;i<n;i++){
		y[i] += h*(c1*k[0][i]+c2*k[1][i]+c3*k[2][i]+c4*k[3][i]+c5*k[4][i]+c6*k[5][i]);
		//printf("\ty[%i] = %g\t",i,y[i]);
		yerr[i] = h*(dc1*k[0][i]+dc2*k[1][i]+dc3*k[2][i]+dc4*k[3][i]+dc5*k[4][i]+dc6*k[5][i]);
		//printf("yerr[%i] = %g\t",i,yerr[i]);
		*maxrelerr = fabs( yerr[i]/(dre*y[i]) );  // (calculated error)/(desired error)   ( Earlier, (dre*y[i])->TOL )
		//printf("*maxrelerr = %g\n",*maxrelerr);
	}

	for(i=0;i<6;i++)  free((char *)k[i]);
	free((char *)ytest);
}



void G_rkf45(int imax, double l[lpts], double dl[lpts], double G[lpts+1], double K[lpts+1], double *t, double *dt, double Gerr[lpts]){
	int i,accurate=0;
	double temprelerr, maxrelerr, factor_dt;
	double *Gtest, *k1,*k2,*k3,*k4,*k5,*k6;  // the k's are slopes, dG/dt

	k1 = (double *)calloc(imax-1,sizeof(double));
	k2 = (double *)calloc(imax-1,sizeof(double));
	k3 = (double *)calloc(imax-1,sizeof(double));
	k4 = (double *)calloc(imax-1,sizeof(double));
	k5 = (double *)calloc(imax-1,sizeof(double));
	k6 = (double *)calloc(imax-1,sizeof(double));
	Gtest = (double *)calloc(imax+1,sizeof(double));
	//Gtemp = (double *)calloc(imax+1,sizeof(double));

	while(accurate==0){
		// get k1 and first Gtest
		for(i=1; i<=imax-1; i++){
			fpGRecRel(l[i], dl[i-1], dl[i], G, K, &k1[i]);
			Gtest[i] = G[i] + *dt*(b21*k1[i]);
		}
		Gtest[0]      = Gtest[1];
		Gtest[imax] = Gtest[imax-1];
		// get k2 and next Gtest
		for(i=1; i<=imax-1; i++){
			fpGRecRel(l[i], dl[i-1], dl[i], Gtest, K, &k2[i]);
			Gtest[i] = G[i] + *dt*(b31*k1[i]+b32*k2[i]);
		}
		Gtest[0]      = Gtest[1];
		Gtest[imax] = Gtest[imax-1];
		// get k3 and next Gtest
		for(i=1; i<=imax-1; i++){
			fpGRecRel(l[i], dl[i-1], dl[i], Gtest, K, &k3[i]);
			Gtest[i] = G[i] + *dt*(b41*k1[i]+b42*k2[i]+b43*k3[i]);
		}
		Gtest[0]      = Gtest[1];
		Gtest[imax] = Gtest[imax-1];
		// get k4 and next Gtest
		for(i=1; i<=imax-1; i++){
			fpGRecRel(l[i], dl[i-1], dl[i], Gtest, K, &k4[i]);
			Gtest[i] = G[i] + *dt*(b51*k1[i]+b52*k2[i]+b53*k3[i]+b54*k4[i]);
		}
		Gtest[0]      = Gtest[1];
		Gtest[imax] = Gtest[imax-1];
		// get k5 and next Gtest
		for(i=1; i<=imax-1; i++){
			fpGRecRel(l[i], dl[i-1], dl[i], Gtest, K, &k5[i]);
			Gtest[i] = G[i] + *dt*(b61*k1[i]+b62*k2[i]+b63*k3[i]+b64*k4[i]+b65*k5[i]);
		}
		Gtest[0]      = Gtest[1];
		Gtest[imax] = Gtest[imax-1];
		// get k6 etc
		for(i=1; i<=imax-1; i++){
			fpGRecRel(l[i], dl[i-1], dl[i], Gtest, K, &k6[i]);
			Gerr[i] = *dt*(dc1*k1[i]+dc2*k2[i]+dc3*k3[i]+dc4*k4[i]+dc5*k5[i]+dc6*k6[i]);
			temprelerr = fabs( Gerr[i]/TOL );  // (calculated error)/(desired error)   ( or TOL->(dre*Gtemp[i]) )
			if(temprelerr>maxrelerr)  maxrelerr = temprelerr;
		}

		if(maxrelerr>1.0){  // calculated error is too big
			accurate = 0;
			// Failure print out //
			printf("G failure: maxrelerr = %g\n",maxrelerr);
			//printf("rkfail: %g\t%g\t%g\t%g\t%g\t%g\t%g\t%i\n", *x, *htemp, ytemp[0], yerr[0], ytemp[1], yerr[1], maxerr, accurate);
			*dt *= SAFETY*pow(maxrelerr,PSHRINK);
		}
		else{  // calculated error small enough (growing next dt)
			accurate = 1;
			// accept dt and calculate next time t and G
			*t += *dt;
			for(i=1; i<=imax-1; i++)  G[i] += *dt*(c1*k1[i]+c2*k2[i]+c3*k3[i]+c4*k4[i]+c5*k5[i]+c6*k6[i]);
			G[0]    = G[1];
			G[imax] = G[imax-1];
			// grow next dt
			factor_dt = SAFETY*pow(maxrelerr,PGROW);
			if(factor_dt<MAXFACTOR)  *dt *= factor_dt;
			else                     *dt *= MAXFACTOR;
		}
	}
	free((char *)k1);
	free((char *)k2);
	free((char *)k3);
	free((char *)k4);
	free((char *)k5);
	free((char *)k6);
	free((char *)Gtest);
	//free((char *)Gtemp);
}



void EqRecRel(double x, double z[2], double dzdx[2], unsigned n){
	dzdx[0] = -B*z[1]*z[1]*z[0]*z[0];
	dzdx[1] = (2-PI*z[0])*z[1];
}



void fpKRecRel(double x, double Kcalc[1], double dKdx[1], unsigned n){
	dKdx[0] = -4*PI*PI*PI*G[i-1]*Kcalc[0]*Kcalc[0]*exp(4*x);
}



void fpGRecRel(double x, double dx1,  double dx2, double Gcalc[lpts+1], double K[lpts+1], double *dGdt){
	// Note that K is the new "t+dt" version of K, and Gcalc is the old "t" version of G
	// Original, constant dx, weird (factors of 2*PI) formula
	//dGdt[0] = exp(-2.0*x)*((1/2.0/PI)*(Gcalc[i+1] - 2.0*Gcalc[i] + Gcalc[i-1])/(dx*dx) + (K[i+1]*Gcalc[i+1] - K[i-1]*Gcalc[i-1])/2.0/dx);
	// Asymmetric derivatives, weird formula
	*dGdt = exp(-2*x)*(  (1/2.0/PI)*( (Gcalc[i+1]-Gcalc[i])/dx2 - (Gcalc[i]-Gcalc[i-1])/dx1 )/dx2  +  (K[i+1]*Gcalc[i+1]-K[i]*Gcalc[i])/dx2  );
	//not dGdt[0] = exp(-2*x)*(  (1/2.0/PI)*( (Gcalc[i+1]-Gcalc[i])/dx2 - (Gcalc[i]-Gcalc[i-1])/dx1 )/dx2  +  (K[i+1]*Gcalc[i+1]-K[i]*Gcalc[i])/dx2  );
	// Asymmetric derivatives, correct formula
	//dGdt[0] = exp(-2*x)*(  2*PI*(K[i+1]*Gcalc[i+1]-K[i]*Gcalc[i])/dx2  +  ( (Gcalc[i+1]-Gcalc[i])/dx2 - (Gcalc[i]-Gcalc[i-1])/dx1 )/dx2  );
	// Symmetric derivatives, weird formula
	//dGdt[0] = exp(-2*x)*(  (1/2.0/PI)*( (Gcalc[i+2]-Gcalc[i])/(dxc+dxd) - (Gcalc[i]-Gcalc[i-2])/(dxa+dxb) )/(dxb+dxc)  +  (K[i+1]*Gcalc[i+1]-K[i-1]*Gcalc[i-1])/(dxb+dxc)  );
	// Symmetric derivatives, correct formula
	//...

	//dzdt = exp(-2*l)*( (zhi[0]*zhi[1]-zlo[0]*zlo[1])/2.0/Dl + (1/2.0/PI)*(zhi[1]-2*zmd[1]+zlo[1])/(Dl*Dl) );
	//dzdt = exp(-2*l)*( PI*(zhi[0]*zhi[1]-zlo[0]*zlo[1])/Dl + (zhi[1]-2*zmd[1]+zlo[1])/(Dl*Dl) );

	// MIGHT HAVE TO FIX THIS SINCE DX+ can be different from DX-
}












/*
		// Before entering loop to get each additional K and G,
		// take l-step to get second K (since each l-step in G requires 3 values of K) //
		my_rkf45(EqRecRel, &l[0], dl[0], Kcalc, 1, DRE, &dltemp, Kcalctemp, Kcalcerr, &maxrelerror);
		dl[0] = dltemp;
		factor = SAFETY*pow(maxrelerror,PGROW);
		if(factor<MAXFACTOR)  dl[1] = dltemp*factor;
		else                  dl[1] = dltemp*MAXFACTOR;
		Kcalc[0] = Kcalctemp[0];
		K[1] = Kcalc[0];
		Kerr[1] = Kcalcerr[0];

SEE vlt_2_ProgramNotes.txt

		get K[i] and dl[i-1]
		get G[i-1] ( with K[i-2],K[i-1],K[i], dl[i-2],dl[i-1] )
			if time-directed error in G[i-1] is too big,
				shrink dt, loop back to calc G[0] to G[i-1] again
				grow next dt
			if  l-directed  error  in G[i-1] is too big,
				shrink dl[i-1]
				if dl[i-1]<0.001dl[i-2]  reduce dl[i-2] by 0.5
				  and [loop back by one]
					recalc K[i-1] and dl[i-1] (and recalc and check G[i-2])
					recalc K[i]
					loop back to calc G[i-1]
			else grow next dt
				record max error

get G[i-1]:     dGdt[i-1] depends on dt and 9 dl[]'s!!!

		// Calculate dGdt for all length scales, and update G //
		fpGRecRel(l[i], dl, G);
		fpGRecRel(get dGdt[i]; x,dx1,dx2, G[i-1],G[i],G[i+1], Ktemp[i],Ktemp[i+1] );
		k1        = dGdt[i]*dt;
		testG1[i] = G[i] + k1/2.0;

		fpGRecRel(l[i], dl, testG1);
		fpGRecRel(get dGdt[i]; x,dx1,dx2, testG1[i-1],testG1[i],testG1[i+1], Ktemp[i],Ktemp[i+1] )
		k2[i]     = dGdt[0]*dt;
		testG2[i] = G[i] - k1[i] + k2[i];

		fpGRecRel(l[i], dl, testG2);
		fpGRecRel(get dGdt[i]; x,dx1,dx2, testG2[i-1],testG2[i],testG2[i+1], Ktemp[i],Ktemp[i+1] )
		k3 = dGdt[0]*dt;
		Gtemp[i]  = G[i] + k1/6.0 + 2.0*k2/3.0 + k3/6.0;




		// Calculate state-variables at progressively larger length scales, up to coherence length, via Runge-Kutta method //
		i = 2; // for stepping l (i = "lstep")
		j = 1; // for selecting data to print to file
		k = 1; // for data labelling
		while(l[i]<lmax){
			// Increment to next lstep, with unknown values for l, K, y, G //
			i++;
////			printf("\ni = %i\n\n",i);
			// Step out in length scale l, calculating next K, y, l, G (and errors in K and G) //
			accuracy = 0;
			while(accuracy==0){
				my_rkf45(EqRecRel, &l[i-1], dl[i-1], Kcalc, 1, DRE, &dltemp, Kcalctemp, Kcalcerr, &maxrelerror);
				// equil: rkf6 2D (K and y), EqRecRel, adaptive step to desired relative error DRE
				//printf("%s\n","made it out of rk!");
				//printf("\tl[%i-1] = %g\n",i,l[i-1]);
				//printf("\tdltemp = %g\n",dltemp);
				l[i] = l[i-1]+dltemp;
				//printf("\tl[%i] = %g\n",i,l[i]);
				//printf("\tztemp[1] = %g\n",ztemp[1]);
				//printf("\texp(-4.0*l[%i]) = %g\n",i,exp(-4.0*l[i]));
				//printf("\tztemp[1]*ztemp[1]*exp(-4.0*l[%i])/a04 = %g\n",i,ztemp[1]*ztemp[1]*exp(-4.0*l[i])/a04);
				G[i] = ztemp[1]*ztemp[1]*exp(-4.0*l[i])/a04;
				//printf("\tG[%i] = %g\n",i,G[i]);
				Gerr[i] = zerr[1]*2*ztemp[1]*exp(-4.0*l[i])/a04;
				//printf("\tGerr[%i] = %g\n",i,Gerr[i]);
				// check error (in G), redo (w/ smaller increment) or proceed (w/o changing increment) //
				Grelerror = fabs( Gerr[i]/TOL );  // (calculated error)/(desired error)   ( Earlier, TOL->(DRE*G[i]) )
				if(Grelerror>maxrelerror)  maxrelerror = Grelerror;
				//printf("\terrfrac = %g\n",errfrac);
				if(maxrelerror>1.0){  // calculated error is too big
////					printf("%s","accuracy = 0\n");
					accuracy = 0;
					dl[i-1] = SAFETY*dltemp*pow(maxrelerror,PSHRINK);
					// Failure print out //
					//printf("G accuracy Fail: %g\t%g\t%g\t%g\n", l[i-1], dl[i-1], ztemp[0], G[i]);
				}
				else{  // calculated error small enough (or overly small), so accept current dl and grow it for next step
					//printf("%s","accuracy = 1\n");
					accuracy = 1;
					dl[i-1] = dltemp;
					factor = SAFETY*pow(maxrelerror,PGROW);
					if(factor<MAXFACTOR)  dl[i] = dltemp*factor;
					else                  dl[i] = dltemp*MAXFACTOR;
					if(l[i]+dl[i]>lmax)   dl[i] = lmax-l[i];
					z[0] = ztemp[0];
					K[i] = z[0];
					Kerr[i] = zerr[0];
					z[1] = ztemp[1];
					// Success print out //
////					printf("G accuracy Success: %g\t%g\t%g\t%g\n", l[i-1], dl[i-1], ztemp[0], G[i]);
				}
			}
			// Print at most [ldatamax] (e.g., 50) data points to screen and output file //
			if(l[i]>=j*dataDl){
				while(l[i]>=j*dataDl)  j++;
				printf(         "%i\t%i\t%g\t%g\t%g\t%g\t%g\t%g\n", k, i, iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, G[i]);
				fprintf(outfile,"%i\t%i\t%g\t%g\t%g\t%g\t%g\t%g\n", k, i, iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, G[i]);
				k++;
			}
		}
		imax = i;

		// Calculate and print total error, min, and max final values //
		for(i=0;i<=imax;i++){
			GerrTot += Gerr[i];
			KerrTot += Kerr[i];
		}
		fprintf(outfile,"%s\t%i\t%g\t%g\t%g\t%g\t%g\t%g\n", "min", imax, iTfrac, 1-iTfrac, l[imax], (K[imax]-KerrTot), (K[imax]-KerrTot)/K0, (G[imax]-GerrTot) );
		fprintf(outfile,"%s\t%i\t%g\t%g\t%g\t%g\t%g\t%g\n", "max", imax, iTfrac, 1-iTfrac, l[imax], (K[imax]+KerrTot), (K[imax]+KerrTot)/K0, (G[imax]+GerrTot) );
		fprintf(outfile,"# Total error in K = %g\n", KerrTot);
		fprintf(outfile,"# Total error in G = %g\n", GerrTot);
		fflush(outfile);

		// Boundary condition enforcement for next stage of program //
		K[imax] = G[imax] = 0.0;  // K=G=0 at l=lmax+dl, terminates recursion relations at l=lmax


	}
*/














//  Program Notes  //
//=========================================================================//
/*

See further documentation (UCLAresearchNotes.pdf, UCLAresearchPrograms.pdf, derk1.c) for elaboration.






int derk(void (*f)(double, double*, double*, unsigned), double *x, double y[], double h, unsigned n){
	int i,j;
	double *k[4],*s,temp;

	for(i=0;i<4;i++)  k[i] = (double *)calloc(n,sizeof(double));
	s = (double *)calloc(n,sizeof(double));
	f(*x,y,k[0],n);
	for(i=1;i<4;i++){
		temp = h*((i+1)/2)/2;
		for(j=0;j<n;j++)  s[j] = y[j]+k[i-1][j]*temp;
		f((*x)+h*(i/2)/2,s,k[i],n);
	}
	for(j=0;j<n;j++)  y[j] += h*(k[0][j]+2*(k[1][j]+k[2][j])+k[3][j])/6;
	*x += h;
	for(i=0;i<4;i++)  free((char *) k[i]);
	free((char *)s);
	return(0);
}

ISN'T THIS INCORRECT?  SHOULDN'T IT BE  f((*x)+h*((i+1)/2)/2,s,k[i],n);


// rk4
const double b21=0.5,
	     b31=0.0,    b32=0.5,
	     b41=0.0,    b42=0.0,    b43=1.0,
	     a1=0.0,     a2=0.5,     a3=0.5,     a4=1.0,
	     c1=1.0/6.0, c2=1.0/3.0, c3=1.0/3.0, c4=1.0/6.0;


void my_rk4(void (*derivs)(double, double*, double*, unsigned), double *x, double y[], double h, unsigned n){
	int i;
	double *k[4],*ytemp;

	for(i=0;i<4;i++)  k[i] = (double *)calloc(n,sizeof(double));
	ytemp = (double *)calloc(n,sizeof(double));

	derivs(*x,y,k[0],n);           // get k1=k[0]
	for(i=0;i<n;i++)  ytemp[i] = y[i]+h*b21*k[0][i];
	derivs(*x+a2*h,ytemp,k[1],n);  // get k2=k[1]
	for(i=0;i<n;i++)  ytemp[i] = y[i]+h*b32*k[1][i];
	derivs(*x+a3*h,ytemp,k[2],n);  // get k3=k[2]
	for(i=0;i<n;i++)  ytemp[i] = y[i]+h*b43*k[2][i];
	derivs(*x+a4*h,ytemp,k[3],n);  // get k4=k[3]

	for(i=0;i<n;i++)  y[i] += h*(c1*k[0][i]+c2*k[1][i]+c3*k[2][i]+c4*k[3][i]);
	*x += h;

	for(i=0;i<4;i++)  free((char *) k[i]);
	free((char *)ytemp);
}


// more general rk4
void my_rk4(void (*derivs)(double, double*, double*, unsigned), double *x, double y[], double h, unsigned n){
	int i;
	double *k[4],*ytemp;

	for(i=0;i<4;i++)  k[i] = (double *)calloc(n,sizeof(double));
	ytemp = (double *)calloc(n,sizeof(double));

	derivs(*x,y,k[0],n);           // get k1=k[0]
	for(i=0;i<n;i++)  ytemp[i] = y[i]+h*b21*k[0][i];
	derivs(*x+a2*h,ytemp,k[1],n);  // get k2=k[1]
	for(i=0;i<n;i++)  ytemp[i] = y[i]+h*(b31*k[0][i]+b32*k[1][i]);
	derivs(*x+a3*h,ytemp,k[2],n);  // get k3=k[2]
	for(i=0;i<n;i++)  ytemp[i] = y[i]+h*(b41*k[0][i]+b42*k[1][i]+b43*k[2][i]);
	derivs(*x+a4*h,ytemp,k[3],n);  // get k4=k[3]

	for(i=0;i<n;i++)  y[i] += h*(c1*k[0][i]+c2*k[1][i]+c3*k[2][i]+c4*k[3][i]);
	*x += h;

	for(i=0;i<4;i++)  free((char *) k[i]);
	free((char *)ytemp);
}


// more general rkf6
void my_rkf6(void (*derivs)(double, double*, double*, unsigned), double *x, double y[], double h, unsigned n){
	int i;
	double *k[6],*ytemp;

	for(i=0;i<6;i++)  k[i] = (double *)calloc(n,sizeof(double));
	ytemp = (double *)calloc(n,sizeof(double));

	derivs(*x,y,k[0],n);           // get k1=k[0]
	for(i=0;i<n;i++)  ytemp[i] = y[i]+h*b21*k[0][i];
	derivs(*x+a2*h,ytemp,k[1],n);  // get k2=k[1]
	for(i=0;i<n;i++)  ytemp[i] = y[i]+h*(b31*k[0][i]+b32*k[1][i]);
	derivs(*x+a3*h,ytemp,k[2],n);  // get k3=k[2]
	for(i=0;i<n;i++)  ytemp[i] = y[i]+h*(b41*k[0][i]+b42*k[1][i]+b43*k[2][i]);
	derivs(*x+a4*h,ytemp,k[3],n);  // get k4=k[3]
	for(i=0;i<n;i++)  ytemp[i] = y[i]+h*(b51*k[0][i]+b52*k[1][i]+b53*k[2][i]+b54*k[3][i]);
	derivs(*x+a5*h,ytemp,k[4],n);  // get k5=k[4]
	for(i=0;i<n;i++)  ytemp[i] = y[i]+h*(b61*k[0][i]+b62*k[1][i]+b63*k[2][i]+b64*k[3][i]+b65*k[4][i]);
	derivs(*x+a6*h,ytemp,k[5],n);  // get k6=k[5]

	for(i=0;i<n;i++)  y[i] += h*(c1*k[0][i]+c2*k[1][i]+c3*k[2][i]+c4*k[3][i]+c5*k[4][i]+c6*k[5][i]);
	*x += h;

	for(i=0;i<6;i++)  free((char *) k[i]);
	free((char *)ytemp);
}


 c1=16.0/135.0,     c2=0.0,             c3=6656.0/12825.0,  c4=28561.0/56430.0,  c5=-9.0/50.0,   c6=2.0/55.0,
 d1=25.0/216.0,     d2=0.0,             d3=1408.0/2565.0,   d4=2197.0/4104.0,    d5=-1.0/5.0,    d6=0.0,
dc1=1.0/360.0,     dc2=0.0,            dc3=-128.0/4275.0,  dc4=-2197.0/75240.0, dc5=1.0/50.0,   dc6=2.0/55.0;
dc1 = (16*216-25*135)/(135*216) = 81/29160 = 1/360
dc2 = 0
dc3 = (6656*2565-1408*12825)/(12825*2565) = -984960/32896125 = -196992/6579225 = -65664/2193075 = -2432/81225 = -128/4275
dc4 = (28561*4104-2197*56430)/(56430*4104) = -6762366/231588720 = -3381183/115794360 = -41743/1429560 = -2197/75240
dc5 = (-9*5+1*50)/(50*5) = 5/(50*5) = 1/50
dc6 = 2/55













//=====================================

y'    = f(t,y)
y(t0) = y0

// kn = test slope = test df/dt
k1 = f(t_n,y_n)
k2 = f(t_n+(1/2)*h,y_n+(1/2)*h*k1)
k3 = f(t_n+(1/2)*h,y_n+(1/2)*h*k2)
k4 = f(t_n+h,y_n+h*k3)

y_{n+1} = y_n + (1/6)*h*(k1+2*k2+2*k3+k4)
t_{n+1} = t_n + h

//=====================================

k1 = f(t_n,y_n)
k2 = f(t_n+c2*h,y_n+a21*h*k1)
k3 = f(t_n+c3*h,y_n+a31*h*k1+a32*h*k2)
...
ks = f(t_n+cs*h,y_n+as1*h*k1+...+a{s,s-1}*h*k{s-1})

y_{n+1} = y_n + h*Sum_{i=1}^s(bi*ki)

Note: Sum_{j=1}^{i-1}aij = ci
      for i = 2,...,s

lower order:
y#_{n+1} = y_n + h*Sum_{i=1}^s(b#i*ki)

er_{n+1} = y_{n+1} - y#_{n+1} = h*Sum_{i=1}^s (bi-b#i)*ki

//=====================================

k1 = f(t_n,y_n)
  ytemp = y_n+a21*h*k1
k2 = f(t_n+c2*h,ytemp)
  ytemp = y_n+a31*h*k1+a32*h*k2
k3 = f(t_n+c3*h,ytemp)
...
  ytemp = y_n+as1*h*k1+...+a{s,s-1}*h*k{s-1}
ks = f(t_n+cs*h,ytemp)

y_{n+1} = y_n + h*Sum_{i=1}^s(bi*ki)

Note: Sum_{j=1}^{i-1}aij = ci
      for i = 2,...,s

lower order:
y#_{n+1} = y_n + h*Sum_{i=1}^s(b#i*ki)

er_{n+1} = y_{n+1} - y#_{n+1} = h*Sum_{i=1}^s (bi-b#i)*ki

//=====================================










a21=0.5;
a31=0.0; a32=0.5;
a41=0.0; a42=0.0; a43=1.0;
c1=0.0; c2=0.5; c3=0.5; c4=1.0;
b1=1.0/6.0; b2=1.0/3.0; b3=1.0/3.0; b4=1.0/6.0;

int derk(void (*derivs)(double, double*, double*, unsigned), double *x, double y[], double h, unsigned n){
	int i,j;
	double *k[4],*s,temp;

	for(i=0;i<4;i++)  k[i] = (double *)calloc(n,sizeof(double));
	s = (double *)calloc(n,sizeof(double));

	derivs(*x,y,k[0],n);

	temp = h*((1+1)/2)/2;				temp = 0.5*h;					
	for(j=0;j<n;j++)  s[j] = y[j]+k[1-1][j]*temp;	for(j=0;j<n;j++)  s[j] = y[j]+k[0][j]*temp;	
	derivs((*x)+h*(1/2)/2,s,k[1],n);		derivs((*x)+0.0*h,s,k[1],n)			

	temp = h*((2+1)/2)/2;				temp = 0.5*h;					
	for(j=0;j<n;j++)  s[j] = y[j]+k[2-1][j]*temp;	for(j=0;j<n;j++)  s[j] = y[j]+k[1][j]*temp;	
	derivs((*x)+h*(2/2)/2,s,k[2],n);		derivs((*x)+0.5*h,s,k[2],n);			

	temp = h*((3+1)/2)/2;				temp = 1.0*h;					
	for(j=0;j<n;j++)  s[j] = y[j]+k[3-1][j]*temp;	for(j=0;j<n;j++)  s[j] = y[j]+k[2][j]*temp;	
	derivs((*x)+h*(3/2)/2,s,k[3],n);		derivs((*x)+0.5*h,s,k[3],n);			

	for(j=0;j<n;j++)  y[j] += h*(k[0][j]+2*(k[1][j]+k[2][j])+k[3][j])/6;
	*x += h;
	for(i=0;i<4;i++)  free((char *) k[i]);
	free((char *)s);
	return(0);
}


a21=0.25;
a31=3.0/32.0;      a32=9.0/32.0;
a41=1932.0/2197.0; a42=-7200.0/2197.0; a43=7296.0/2197.0;
a51=439.0/216.0;   a52=-8.0;           a53=3680.0/513.0;   a54=-845.0/4104.0;
a61=-8.0/27.0;     a62=2.0;            a63=-3544.0/2565.0; a64=1859.0/4104.0; a65=-11.0/40.0;
c1=0.0;        c2=0.25; c3=0.375;          c4=12.0/13.0;       c5=1.0;        c6=0.5;
b1=16.0/135.0; b2=0.0;  b3=6656.0/12825.0; b4=28561.0/56430.0; b5=-9.0/50.0;  b6=2.0/55.0;
d1=25.0/216.0; d2=0.0;  d3=1408.0/2565.0;  d4=2197.0/4104.0;   d5=-1.0/5.0;   d6=0.0;


k1 = f(t_n,y_n)
  ytemp = y_n+a21*h*k1
k2 = f(t_n+c2*h,ytemp)
  ytemp = y_n+a31*h*k1+a32*h*k2
k3 = f(t_n+c3*h,ytemp)
  ytemp = y_n+a41*h*k1+a42*h*k2+a43*h*k3
k4 = f(t_n+c4*h,ytemp)
  ytemp = y_n+a51*h*k1+a52*h*k2+a53*h*k3+a54*h*k4
k5 = f(t_n+c5*h,ytemp)
  ytemp = y_n+a61*h*k1+a62*h*k2+a63*h*k3+a64*h*k4+a65*h*k5
k6 = f(t_n+c6*h,ytemp)

y_{n+1}  = y_n + h*Sum_{i=1}^s(bi*ki)
y#_{n+1} = y_n + h*Sum_{i=1}^s(b#i*ki)

	ytemp=vector(1,n);
	for(i=1;i<=n;i++)  // First step.
		ytemp[i]=y[i]+b21*h*dydx[i];
	(*derivs)(x+a2*h,ytemp,ak2);  // Second step.
	for(i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
	(*derivs)(x+a3*h,ytemp,ak3);  // Third step.
	for(i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
	(*derivs)(x+a4*h,ytemp,ak4);  // Fourth step.
	for(i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
	(*derivs)(x+a5*h,ytemp,ak5);  // Fifth step.
	for(i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
	(*derivs)(x+a6*h,ytemp,ak6);  // Sixth step.
	for(i=1;i<=n;i++)  // Accumulate increments with proper weights.
		yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
	for(i=1;i<=n;i++)
		yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
		// Estimate error as difference between fourth and fifth order methods.




// more general rkf6
void my_rkf6(){
	int i;
	double *k[6],*s;

	for(i=0;i<6;i++)  k[i] = (double *)calloc(n,sizeof(double));
	s = (double *)calloc(n,sizeof(double));

	derivs(*x,y,k1,n);           // get k1
	ytemp = y+b21*h*k1;
	derivs(*x+a2*h,ytemp,k2,n);  // get k2
	ytemp = y+b31*h*k1+b32*h*k2;
	derivs(*x+a3*h,ytemp,k3,n);  // get k3
	ytemp = y+b41*h*k1+b42*h*k2+b43*h*k3;
	derivs(*x+a4*h,ytemp,k4,n);  // get k4
	ytemp = y+b51*h*k1+b52*h*k2+b53*h*k3+b54*h*k4;
	derivs(*x+a5*h,ytemp,k5,n);  // get k5
	ytemp = y+b61*h*k1+b62*h*k2+b63*h*k3+b64*h*k4+b65*h*k5;
	derivs(*x+a6*h,ytemp,k6,n);  // get k6

	for(i=0;i<n;i++)  y[i] += h*(c1*k[0][i]+c2*k[1][i]+c3*k[2][i]+c4*k[3][i]+c5*k[4][i]+c6*k[5][i]);
	*x += h;

	for(i=0;i<6;i++)  free((char *) k[i]);
	free((char *)s);
}




*/
