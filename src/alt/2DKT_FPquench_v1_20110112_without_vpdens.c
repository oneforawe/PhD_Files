//  File Comments  //
//=========================================================================//

/* FILENAME: 2DKT_FPquench.c */
/* VERSION: 1 (2011 Jan 12)
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

// Program parameters/inputs and data type definitions //
const int    lmax     = 10;        //  1  //  maximum of the length scale l (potentially approximate)
const int    lsteps   = 1250;      //1250 // 5000 // 100000 //  lsteps = lmax/dl = 10/0.0001 = 100,000  (make sure this is divisible by ldatamax)
const int    lpts     = lsteps+1;  //  from l=0 to l=lmax, inclusive (incl. the boundary-condition-enforcing pnt, there are lpts+1 "l" pnts)
const int    ldatamax = 50;        //  100 //  max number of data points recorded per time step
const double iTfrac   = 0.95;      //  initial temperature fraction Ti/Tkt
const double qTfrac   = 0.10;      //  quench temperature fraction Tq/Tkt
const double tmax     = 10000.0;   //  20600 //  max unitless time (see if(t==...) below)
const double dt0      = 1.0e-5;    //  the time increment (in units of the "diffusion time")
const double a0       = 1.0;       //  a0 in units of a0 (!)
const double a04      = 1.0;       //  a0 to the fourth power
const double K0c      = 0.747853;  //  critical (T=Tkt) bare superfluid ratio (or "coupling constant")

int i;
double K[lpts+1], G[lpts+1];
double dGdt[1];



//  Function prototypes  //
//=========================================================================//
// derk((*f)(),*x,y[],h,n);
// EqRecRel(x,z,dzdx,n);
// fpKRecRel(x,Kcalc,dKdx,n);
// fpGRecRel(x,Dl,Gcalc);
int  derk(void (*f)(double, double*, double*, unsigned int), double *x, double y[], double h, unsigned n);
void EqRecRel(double x, double z[2], double dzdx[2], unsigned n);
void fpKRecRel(double x, double Kcalc[1], double dKdx[1], unsigned n);
void fpGRecRel(double x, double Dl, double Gcalc[lsteps+1]);



//  Function definitions  //
//=========================================================================//

main(){
	// Main function definitions //
	int64_t n;
	int k=1;
	double dblmax=lmax, dblsteps=lsteps;
	double K0, y0,  dl=dblmax/dblsteps;
	double x, l[lpts], z[2];
	double Kcalc[lpts+1], Gcalc[lpts+1];
	double k1[lpts], k2[lpts], k3[lpts];
	double testG1[lpts+1], testG2[lpts+1];
	double t, dt, dtLimit;
	FILE *outfile;
	char *filename;

	// Prepare output file, print identification and values //
	asprintf(&filename, "2DKT_FPquench_lmax%i_dl%g_tmax%g_dt0%g.dat", lmax, dl, tmax, dt0);
	outfile = fopen(filename,"w");  //  E.g., "2DKT_FPquench_adapt.out"
	fprintf(outfile,"# Filename: %s\n", filename);
	fprintf(outfile,"# Source: 2DKT_FPquench.c\n");
	fprintf(outfile,"# Source version: %s\n", "1 (2010 12 09)");
	fprintf(outfile,"# Parameter values: lmax=%i, lsteps=%i, lpts=%i, ldatamax=%i, iTfrac=%g, qTfrac=%g, tmax=%g, dt0=%g, a0=%g, a04=%g, K0c=%g\n", lmax,lsteps,lpts,ldatamax,iTfrac,qTfrac,tmax,dt0,a0,a04,K0c);

	// Boundary condition enforcement //
	K[lpts] = G[lpts] = 0;  // K=G=0 at l=lmax+dl, terminates recursion relations at l=lmax


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
	x = 0.0;
	l[0] = 0.0;
	K0 = K0c/iTfrac;
	y0 = exp(-PISQ*K0/2);
	z[0] = K0;
	z[1] = y0;
	K[0] = z[0];
	G[0] = z[1]*z[1]*exp(-4.0*l[0])/a04;  // G (not = 0.0)

	// Print initial data to screen and output file //
	printf(         "\n# time step n = %ld\tt = %e\t(dt = %e)\n", n,t,dt);
	fprintf(outfile,"\n# time step n = %ld\tt = %e\t(dt = %e)\n", n,t,dt);
	printf(         "%s\t%s\t%s\t%s\t%s\t%s\n",   "T/Tc","1-T/Tc","l","K","K/K0","G");
	fprintf(outfile,"# %s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","G");
	printf(         "%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[0], K[0], K[0]/K0, G[0]);
	fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[0], K[0], K[0]/K0, G[0]);

	// Calculate state-variables at progressively larger length scales, up to coherence length, via Runge-Kutta method (derk) //
	for(i=1; i<lpts; i++){
		// Step out in length scale l, calculating next K, y, l, and G //
		derk(EqRecRel, &x, z, dl, 2);  // equil: rk4 2D (K and y), EqRecRel
		l[i] = l[i-1] + dl;
		K[i] = z[0]; 
		G[i] = z[1]*z[1]*exp(-4.0*l[i])/a04;
		// check errors (in K and y), redo (w/ smaller increment) or proceed (setting larger increment) //
		// ...code
		// Print at most [ldatamax] (e.g., 50) data points to screen and output file //
		if(i%(lsteps/ldatamax)==0){
			printf(         "%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, G[i]);
			fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, G[i]);
		}
	}
	fflush(outfile);


	//////////////////////////////////////////////////////////////
	//                                                          //
	//  (n>=1) After Instantaneous "Quench" (Temperature Drop)  //
	//                                                          //
	//////////////////////////////////////////////////////////////

	// Initialize //
	K0 = K0c/qTfrac;
	//y0 = exp(-PISQ*K0/2.0);
	K[0] = K0;
	//G[0] = exp(-4.0*l[0])*y0*y0/a04;  // G (not = 0.0)  //  (correct location?)

	// Calculate K at progressively larger length scales, up to coherence length, via Runge-Kutta method (derkK) //
	x = 0.0;
	Kcalc[0] = K[0];
	for(i=1; i<lpts; i++){
		derk(fpKRecRel, &x, Kcalc, dl, 1);  // nequil: rk4 1D (K), fpKRecRel
		K[i] = Kcalc[0];
		// check error (in K), redo (w/ smaller increment) or proceed (setting larger increment) //
	}

	// Advance in time: rk3 1D (G), fpGRecRel //
	while(t<=tmax+dt){
		n++;
		t += dt;
		printf(         "# time step n = %ld\tt = %e\t(dt = %e)\t", n,t,dt);
		//fprintf(outfile,"# time step n = %ld\tt = %e\t(dt = %e)\t", n,t,dt);

		// Calculate dGdt for all length scales, and update G //
			// nequil: rk4 1D (K), fpKRecRel //
		for(i=1; i<lpts-1; i++){
			fpGRecRel(l[i], dl, G);
			k1[i]     = dGdt[0]*dt;
			testG1[i] = G[i] + k1[i]/2.0;
		}
		testG1[0]      = testG1[1];
		testG1[lpts-1] = testG1[lpts-2];

		for(i=1; i<lpts-1; i++){
			fpGRecRel(l[i], dl, testG1);
			k2[i]     = dGdt[0]*dt;
			testG2[i] = G[i] - k1[i] + k2[i];
		}
		testG2[0]      = testG2[1];
		testG2[lpts-1] = testG2[lpts-2];

		for(i=1; i<lpts-1; i++){
			fpGRecRel(l[i], dl, testG2);
			k3[i] = dGdt[0]*dt;
			G[i]  = G[i] + k1[i]/6.0 + 2.0*k2[i]/3.0 + k3[i]/6.0;  //  update G here and below**
			// check error (in G), redo (w/ smaller increment) or proceed (setting larger increment) //
		}
		G[0]      = G[1];       //  **here
		G[lpts-1] = G[lpts-2];  //  **and here

		// Initialize //
		K0 = K0c/qTfrac;
		K[0] = K0;  //  K

		// Calculate K at progressively larger length scales, up to coherence length, via Runge-Kutta method (derk1) //
		x = 0.0;
		Kcalc[0] = K[0];
		for(i=1; i<lpts; i++){
			derk(fpKRecRel, &x, Kcalc, dl, 1);  // nequil: rk4 1D (K), fpKRecRel
			K[i] = Kcalc[0];
		}

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

	fclose(outfile);
}



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



void EqRecRel(double x, double z[2], double dzdx[2], unsigned n){
	dzdx[0] = -B*z[1]*z[1]*z[0]*z[0];
	dzdx[1] = (2-PI*z[0])*z[1];
}



void fpKRecRel(double x, double Kcalc[1], double dKdx[1], unsigned n){
	dKdx[0] = -4*PI*PI*PI*G[i-1]*Kcalc[0]*Kcalc[0]*exp(4.0*x);
}



void fpGRecRel(double x, double Dl, double Gcalc[lsteps+1]){
	//dGdt[0] = exp(-2.0*x)*((1/2.0/PI)*(Gcalc[i+1] - 2.0*Gcalc[i] + Gcalc[i-1])/(Dl*Dl) + (K[i+1]*Gcalc[i+1] - K[i-1]*Gcalc[i-1])/2.0/Dl);
	dGdt[0] = exp(-2*x)*( PI*(K[i+1]*Gcalc[i+1]-K[i-1]*Gcalc[i-1])/Dl + (Gcalc[i+1]-2*Gcalc[i]+Gcalc[i-1])/(Dl*Dl) );
	//dzdt = exp(-2*l)*( (zhi[0]*zhi[1]-zlo[0]*zlo[1])/2.0/Dl + (1/2.0/PI)*(zhi[1]-2*zmd[1]+zlo[1])/(Dl*Dl) );
	//dzdt = exp(-2*l)*( PI*(zhi[0]*zhi[1]-zlo[0]*zlo[1])/Dl + (zhi[1]-2*zmd[1]+zlo[1])/(Dl*Dl) );
}



//  Program Notes  //
//=========================================================================//
/*

See further documentation (UCLAresearchNotes.pdf, UCLAresearchPrograms.pdf, derk1.c) for elaboration.




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

lower order

y#_{n+1} = y_n + h*Sum_{i=1}^s(b#i*ki)

er_{n+1} = y_{n+1} - y#_{n+1} = h*Sum_{i=1}^s (bi-b#i)*ki

//=====================================

*/


