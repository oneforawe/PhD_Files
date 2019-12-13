//  File Comments  //
//=========================================================================//

/* FILENAME: 2DKT_FPquench.c */
/* VERSION: 5 (2011 Jun 10 - ...)
            Fixing the fact that Gamma blows up when the step size dl is small (10,000 steps with lmax=10, blows up on fourth time-step).
            Got rid of extra point (lpts+1).
            With multiple options for boundary conditions of Gamma.
            Eqns in terms of Gamma (not y), Gamma should no longer fall below its new equilibrium curve */
/* AUTHOR: Andrew Forrester <aforrester@ucla.edu> */
/* DESCRIPTION:
   * This C program calculates time-dependent properties of a 2D system of superfluid (a thin layer of liquid 4He) before and after an instantaneous quench (temperature drop) starting at temperature Ti (Ti/Tkt = iTfrac) at time t=0 and going to temperature Tq (Tq/Tkt = qTfrac) at the next time step.  The system starts in equilibrium and the quench puts the system out of equilibrium so the properties evolve toward a new equilibrium.
   * The properties calculated are superfluid ratios (K, K/K0 = rho/rho0) and vortex-pair probability density G at various length-scales (l).
*/
/* EXT FILES: none */
/* COMPILE NOTES:
   * To compile, type "g++ -lm 2DKT_FPquench.c" without the quotes; then to run, type "./a.out".
   * NOTE: Be sure lsteps is divisible by ldatamax.
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
const int    BC       = 3;            //  Boundary Condition choice 1=Hanching's, 2=Gary's, 3=Mine
const double iTfrac   = 1.00;         //  initial temperature fraction Ti/Tkt
const double qTfrac   = 0.10;         //  quench temperature fraction Tq/Tkt
const double lmax     = 10.0;         //  1  //  maximum of the length scale l (potentially approximate)
const int    lsteps   = 80000;         //  1250, 5000, 100000 //  lsteps = lmax/dl = 10/0.0001 = 100,000  (make sure this is divisible by ldatamax)
const int    lpts     = lsteps+1;     //  number of "l" points, from l=0 to l=lmax, inclusive
const int    ldatamax = 100;          //  50,100 //  max number of data points recorded per time step
const double t0       = 1.0;          //  the diffusion time ( q/a02*k^2 = Lambda*kB*T/a0^2*kappa^2*sigmas0^2 ) in units of t0 (!)
const double tmax     = 10.0;         //  10000.0 //  max time (in units of the diffusion time)
const double dt0      = 1.0e-5;       //  the time increment (in units of the diffusion time)
const double a0       = 1.0;          //  a0 in units of a0 (!)
const double a02      = a0*a0;        //  a0 to the second power
const double a04      = a02*a02;      //  a0 to the fourth power
const double K0c      = 0.747853;     //  critical (T=Tkt) bare superfluid ratio (or "coupling constant")
/*
// Parameters for Runge-Kutta method, 3rd order (Version 1) RK3v1
const double B21=0.5,
	     B31=-1.0,    B32=1.0,
	      A1=0.0,      A2=0.0,      A3=0.0,
	      C1=1.0/6.0,  C2=2.0/3.0,  C3=1.0/6.0;

// Parameters for Runge-Kutta method, 3rd order (Version 2) RK3v2
const double B21=0.5,
	     B31=-1.0,    B32=2.0,
	      A1=0.0,      A2=0.5,      A3=1.0,
	      C1=1.0/6.0,  C2=2.0/3.0,  C3=1.0/6.0;
*/
// Parameters for Runge-Kutta method, 3rd order (Version 3) RK3v3
const double B21=1.0/3.0,
	     B31=0.0,     B32=2.0/3.0,
	      A1=0.0,      A2=1.0/3.0,  A3=2.0/3.0,
	      C1=1.0/4.0,  C2=0.0,      C3=3.0/4.0;

// Global variables //
int i;
double K[lpts], G[lpts];
double dGdt[1];



//  Function prototypes  //
//=========================================================================//

// rk4((*f)(),*x,y[],h,n);
// EqRecRel(x,z,dzdx,n);
// fpKRecRel(x,Kcalc,dKdx,n);
// fpGRecRel(x,Dl,Gcalc);
void rk4(void (*f)(double, double*, double*, unsigned int), double *x, double y[], double h, unsigned n);
void EqRecRel(double x, double z[2], double dzdx[2], unsigned n);
void fpKRecRel(double x, double Kcalc[1], double dKdx[1], unsigned n);
void fpGRecRel(double x, double Dl, double Gcalc[lpts]);
void fpGRecRelprint(double x, double Dl, double Gcalc[lpts]);



//  Function definitions  //
//=========================================================================//

main(){
	// Main function definitions //
	int64_t n_t;
	int j,k;
	double dblsteps=lsteps, dl=lmax/dblsteps;
	double x, l[lpts];
	double K0, G0, G0final;
	double z[2], Kcalc[lpts], Gcalc[lpts];
	double k1[lpts], k2[lpts], k3[lpts];
	double testG1[lpts], testG2[lpts];
	double t, dt, ttrig;
	double vpintegral, vpdens, vpdens_old;
	double prob_out, prob_dec;
	FILE *outfile1;  // recording functions versus length-scale l ("vsl"), at selected times
	FILE *outfile2;  // recording functions versus time t ("vst"), at largest length scale lmax
	char *filename1, *filename2;

	// Prepare output files, print identification and values //
	asprintf(&filename1, "2DKT_FPquench_vsl_T_%g_%g_lmax%g_dl%g_tmax%g_dt0%g_BC%i_dbug_NewRK1v3.dat", iTfrac, qTfrac, lmax, dl, tmax, dt0, BC);
	asprintf(&filename2, "2DKT_FPquench_vst_T_%g_%g_lmax%g_dl%g_tmax%g_dt0%g_BC%i_dbug_NewRK1v3.dat", iTfrac, qTfrac, lmax, dl, tmax, dt0, BC);
	outfile1 = fopen(filename1,"w");  //  E.g., "2DKT_FPquench.out"
	outfile2 = fopen(filename2,"w");  //  E.g., "2DKT_FPquench.out"
	fprintf(outfile1,"# Filename: %s\n", filename1);
	fprintf(outfile2,"# Filename: %s\n", filename2);
	fprintf(outfile1,"# Source: 2DKT_FPquench.c\n");
	fprintf(outfile2,"# Source: 2DKT_FPquench.c\n");
	fprintf(outfile1,"# Source version: %s\n", "4 (2011 May 26 - ...)");
	fprintf(outfile2,"# Source version: %s\n", "4 (2011 May 26 - ...)");
	fprintf(outfile1,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g, B=%21.21g\n", PI,PISQ,PICU,B);
	fprintf(outfile2,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g, B=%21.21g\n", PI,PISQ,PICU,B);
	fprintf(outfile1,"# Parameter values: BC=%i, iTfrac=%g, qTfrac=%g, lmax=%g, lsteps=%i, (dl=%g), lpts=%i, ldatamax=%i, tmax=%g, dt0=%g, a0=%g, a02=%g, a04=%g, K0c=%g\n", BC,iTfrac,qTfrac,lmax,lsteps,dl,lpts,ldatamax,tmax,dt0,a0,a02,a04,K0c);
	fprintf(outfile2,"# Parameter values: BC=%i, iTfrac=%g, qTfrac=%g, lmax=%g, lsteps=%i, (dl=%g), lpts=%i, ldatamax=%i, tmax=%g, dt0=%g, a0=%g, a02=%g, a04=%g, K0c=%g\n", BC,iTfrac,qTfrac,lmax,lsteps,dl,lpts,ldatamax,tmax,dt0,a0,a02,a04,K0c);

	// Boundary condition enforcement //
	K[lpts] = G[lpts] = 0.0;  // K=G=0 at l=lmax+dl, terminates recursion relations at l=lmax


	////////////////////////////////////////////////
	//                                            //
	//  (n_t=0) Initial Equilibrium Calculations  //
	//                                            //
	////////////////////////////////////////////////

	// Initialize time quantities //
	n_t = 0;
	t = 0.0;
	dt = dt0;  // We won't step in time until the quench.

	// Initialize quantities at smallest length-scale (l=0) //
	x = l[0] = 0.0;
	z[0] = K[0] = K0 = K0c/iTfrac;
	z[1] = G[0] = G0 = exp(-4.0*l[0]-PISQ*K0)/a04;
	vpintegral = 0.0;
	vpintegral += 0.5*G[0]*exp(2.0*l[0])*dl;

	// Print initial data to screen and output files //
	printf(          "\n# time step n_t = %ld\tt = %e\t(dt = %e)\n", n_t,t,dt);
	fprintf(outfile1,"\n# time step n_t = %ld\tt = %e\t(dt = %e)\n", n_t,t,dt);
	printf(          "%s\t%s\t%s\t%s\t%s\t%s\n",   "T/Tc","1-T/Tc","l","K","K/K0","G");
	fprintf(outfile1,"# %s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","G");
	fprintf(outfile2,"# %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "time-step", "t", "T/Tc","1-T/Tc","l","K","K/K0","G","vpdens");
	printf(          "%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[0], K[0], K[0]/K0, G[0]);
	fprintf(outfile1,"%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[0], K[0], K[0]/K0, G[0]);

	// Calculate state-variables at progressively larger length scales, up to coherence length, via Runge-Kutta method (rk4) //
	for(i=1; i<lpts; i++){
		// Step out in length scale l, calculating next l, K, and G //
		rk4(EqRecRel, &x, z, dl, 2);  // equil: rk4 2D (K and G), EqRecRel
		l[i] = l[i-1] + dl;
		K[i] = z[0];
		G[i] = z[1];
		vpintegral += G[i]*exp(2.0*l[i])*dl;
		// Print at most [ldatamax] (e.g., 50) data points to screen and output file 1 //
		if(i%(lsteps/ldatamax)==0 || i<15){
			printf(          "%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, G[i]);
			fprintf(outfile1,"%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, G[i]);
		}
	}
	vpdens = 2.0*PI*a02*vpintegral;
	fflush(outfile1);

	// Print out the data at largest length scale to file 2 //
	fprintf(outfile2,"%ld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", n_t, t, iTfrac, 1-iTfrac, l[lpts-1], K[lpts-1], K[lpts-1]/K0, G[lpts-1], vpdens);
	fflush(outfile2);


	////////////////////////////////////////////////////////////////
	//                                                            //
	//  (n_t>=1) After Instantaneous "Quench" (Temperature Drop)  //
	//                                                            //
	////////////////////////////////////////////////////////////////

	// Initialize //
	K[0] = K0 = K0c/qTfrac;

	// Boundary condition: don't let G0 go below this equilibrium value //
	G0final = exp(-4.0*l[0]-PISQ*K0)/a04;  // = exp(-PISQ*K0)/a04;

	// Calculate K at progressively larger length scales, up to coherence length, via Runge-Kutta method (rk4 K) //
	x = 0.0;
	Kcalc[0] = K[0];
	for(i=1; i<lpts; i++){
		rk4(fpKRecRel, &x, Kcalc, dl, 1);  // nequil: rk4 1D (K), fpKRecRel
		K[i] = Kcalc[0];
	}

	// Advance in time: rk3 1D (G), fpGRecRel //
	j = 1;        // for selecting times to print out l-dependent data: outfile1
	k = 1;        // for selecting times to print out t-dependent data (at largest length scale): outfile2
	ttrig = 0.0;  // for selecting times to print out t-dependent data (at largest length scale): outfile2
	//while(t<=tmax+dt){
	while(n_t<=3){
		n_t++;
		t += dt;
		vpdens_old = vpdens;
		vpintegral = 0.0;
		printf(          "# time step n_t = %ld\tt = %e\t(dt = %e)\t", n_t,t,dt);
		//fprintf(outfile1,"# time step n_t = %ld\tt = %e\t(dt = %e)\t", n_t,t,dt);

		// Calculate dGdt for all length scales, and update G using third order Runge-Kutta method //
			// nequil: rk3 1D (G), fpGRecRel //
		for(i=1; i<lpts-1; i++){
			if(i==1) fpGRecRelprint(l[i]+A1*dl, dl, G);
			else     fpGRecRel(l[i]+A1*dl, dl, G);
			k1[i]     = dGdt[0]*dt;
			testG1[i] = G[i] + B21*k1[i];
		}
		k1[0]     = k1[1];                             //*could be questioned
		testG1[0] = G[0] + B21*k1[0];                  //*
		k1[lpts-1]     = 2.0*k1[lpts-2] - k1[lpts-3];  //*
		testG1[lpts-1] = G[lpts-1] + B21*k1[lpts-1];   //*

		for(i=1; i<lpts-1; i++){
			if(i==1) fpGRecRelprint(l[i]+A2*dl, dl, testG1);
			else     fpGRecRel(l[i]+A2*dl, dl, testG1);
			k2[i]     = dGdt[0]*dt;
			testG2[i] = G[i] + B31*k1[i] + B32*k2[i];
		}
		k2[0]     = k2[1];                                             //*
		testG2[0] = G[0] + B31*k1[0] + B32*k2[0];                      //*
		k2[lpts-1]     = 2.0*k2[lpts-2] - k2[lpts-3];                  //*
		testG2[lpts-1] = G[lpts-1] + B31*k1[lpts-1] + B32*k2[lpts-1];  //*

		for(i=1; i<lpts-1; i++){
			if(i==1) fpGRecRelprint(l[i]+A3*dl, dl, testG2);
			else     fpGRecRel(l[i]+A3*dl, dl, testG2);
			k3[i] = dGdt[0]*dt;
			G[i]  = G[i] + C1*k1[i] + C2*k2[i] + C3*k3[i];
			vpintegral += G[i]*exp(2*l[i])*dl;
		}
		// Large-scale boundary:
		// assuming that the large-scale behavior is not yet affected by the quench,
		// we can use the equilibrium relation  dG/dl = -2piKG  to derive:
		G[lpts-1] = G[lpts-2] - 2.0*PI*K[lpts-2]*G[lpts-2]*dl;
		// or G[lpts-1] = G[lpts-2]/(1.0+2.0*PI*K[lpts-1]*dl);
		vpintegral += G[lpts-1]*exp(2.0*l[lpts-1])*dl;
		vpdens = 2.0*PI*a02*vpintegral;  // temporary incomplete calculation of vpdens
		// Small-scale boundary:
		switch(BC){
			case 1: // Hanching's bndry cond. Second l-derivative of G is set to zero (linear continuation).
				G[0] = 2.0*G[1] - G[2];
				if(G[0]<=G0final)  G[0] = G0final;
				break;
			case 2: // Gary's bndry cond. Immediate and permanent drop to new equilibrium value.
				G[0] = G0final;
				break;
			case 3: // My bndry cond. Total probability outflow must occur at smallest scale.
				G[0] = (t0*(vpdens_old-vpdens)*dl/(2.0*PI*a02) - G[1]*dt) / (2.0*PI*K[0]*dl*dt + 0.5*t0*dl*dl - dt);
				if(G[0]<=G0final)  G[0] = G0final;
				break;
		}
		vpintegral += 0.5*G[0]*exp(2*l[0])*dl;
		vpdens = 2.0*PI*a02*vpintegral;

		// Check to make sure probability outflow matches probability decrease
		prob_out = (a0/t0)*(2.0*PI*G[0]*K[0]+(G[1]-G[0])/dl)*2.0*PI*a0;
		prob_dec = -(vpdens-vpdens_old)/dt;
		//printf("\n\nprob_out = %21.21g\tprob_dec = %21.21g", prob_out, prob_dec);
		if(prob_out-prob_dec > 1.0e-6){
			printf("\n\nprob_out - prob_dec = %g", prob_out-prob_dec);
			printf("\nCode Problem: Vortex-pair probability outflow doesn't equal decrease in probability!\n");
			exit(EXIT_FAILURE);
		}

		// Calculate K at progressively larger length scales, up to coherence length, via Runge-Kutta method (rk4 K) //
		x = 0.0;
		Kcalc[0] = K[0];
		for(i=1; i<lpts; i++){
			rk4(fpKRecRel, &x, Kcalc, dl, 1);  // nequil: rk4 1D (K), fpKRecRel
			K[i] = Kcalc[0];
		}

		// At certain times, print out the results to file 1 //
		//if( (t>=0.01&&j==1) || (t>=0.1&&j==2) || (t>=1&&j==3) || (t>=10&&j==4) || (t>=100&&j==5) || (t>=1000&&j==6) || (t>=10000&&j==7) ){
		if( (n_t==1&&j==1) || (n_t==2&&j==2) || (n_t==3&&j==3) || (n_t==4&&j==4) ){
		//if( (n_t==10&&j==1) || (n_t==20&&j==2) || (n_t==30&&j==3) || (n_t==40&&j==4) || (t>=0.01&&j==5) || (t>=0.1&&j==6) || (t>=1&&j==7) || (t>=10&&j==8) || (t>=100&&j==9) || (t>=1000&&j==10) || (t>=10000&&j==11) ){
		//if(  ){
			fprintf(outfile1,"\n# time step n_t = %ld\tt = %e\t(dt = %e)\n", n_t,t,dt);
			fprintf(outfile1,"# %s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","G");
			for(i=0; i<lpts; i++){
				// Print at most [ldatamax] (e.g., 1000) data points to screen and output file //
				if(i%(lsteps/ldatamax)==0 || i<15){
					printf(          "%g\t%g\t%g\t%g\t%g\t%g\n", qTfrac, 1-qTfrac, l[i], K[i], K[i]/K0, G[i]);
					fprintf(outfile1,"%g\t%g\t%g\t%g\t%g\t%g\n", qTfrac, 1-qTfrac, l[i], K[i], K[i]/K0, G[i]);
				}
			}
			fflush(outfile1);
			j++;
		}

		// Print out the data at the largest length-scale to file 2, selecting data at exponentially sparser time intervals //
		if(t>ttrig){
			fprintf(outfile2,"%ld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", n_t, t, qTfrac, 1-qTfrac, l[lpts-1], K[lpts-1], K[lpts-1]/K0, G[lpts-1], vpdens);
			fflush(outfile2);
			ttrig += dt0*pow(10.0,0.05*k);
			//ttrig += dt0*exp(0.05*k);
			k++;
		}

		// Bug detection //
		printf(          "%g\t%g\t%g\n", l[lpts-1], K[lpts-1]/K0, G[lpts-1]);
		//fprintf(outfile1,"%g\t%g\t%g\n", l[lpts-1], KG[lpts-1].c[0]/K0, KG[lpts-1].c[1]);
	}

	fclose(outfile1);
	fclose(outfile2);
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



// equilibrium (K,G) recursion relations //
void EqRecRel(double x, double z[2], double dzdx[2], unsigned n){
	dzdx[0] = -B*a04*exp(4.0*x)*z[1]*z[0]*z[0];
	dzdx[1] = -2.0*PI*z[0]*z[1];
}



// Fokker-Planck K recursion relations //
void fpKRecRel(double x, double Kcalc[1], double dKdx[1], unsigned n){
	dKdx[0] = -B*a04*exp(4*x)*G[i-1]*Kcalc[0]*Kcalc[0];
}



// Fokker-Planck G recursion relations //
void fpGRecRel(double x, double Dl, double Gcalc[lpts]){
	// My version of the equation:
	dGdt[0] = exp(-2*x)*( PI*(K[i+1]*Gcalc[i+1]-K[i-1]*Gcalc[i-1])/Dl + (Gcalc[i+1]-2*Gcalc[i]+Gcalc[i-1])/(Dl*Dl) );
	// Hanching Chu's version, which is mine divided by 2*PI (AND k is multiplied with the unmodified g, as opposed to gin):
	//dGdt[0] = exp(-2.0*x)*( (K[i+1]*Gcalc[i+1]-K[i-1]*Gcalc[i-1])/2.0/Dl + (1/2.0/PI)*(Gcalc[i+1]-2.0*Gcalc[i]+Gcalc[i-1])/(Dl*Dl) );
	//gstep[0]= exp(-2.0*x)*( (k[j+1]*g[j+1]-k[j-1]*g[j-1])/2.0/stepsize + (1/2.0/PI)*(gin[j+1]-2.0*gin[j]+gin[j-1])/(stepsize*stepsize) );
}


/*
// Fokker-Planck G recursion relations //
void fpGRecRelmidstep(double x, double Dl, double Gcalc[lpts]){
	dGdt[0] = exp(-2*x)*( PI*(K[i+1]*Gcalc[i+1]-K[i-1]*Gcalc[i-1])/Dl + (Gcalc[i+1]-2*Gcalc[i]+Gcalc[i-1])/(Dl*Dl) );
}
*/


// Fokker-Planck G recursion relations //
void fpGRecRelprint(double x, double Dl, double Gcalc[lpts]){
	dGdt[0] = exp(-2*x)*( PI*(K[i+1]*Gcalc[i+1]-K[i-1]*Gcalc[i-1])/Dl + (Gcalc[i+1]-2*Gcalc[i]+Gcalc[i-1])/(Dl*Dl) );
	printf("dGdt[0]\n = exp(-2*x)*( PI*(K[i+1]*Gcalc[i+1]-K[i-1]*Gcalc[i-1])/Dl + (Gcalc[i+1]-2*Gcalc[i]+Gcalc[i-1])/(Dl*Dl) )\n = (%g)*( (%g)*((%g)*(%g)-(%g)*(%g))/(%g) + ((%g)-2*(%g)+(%g))/(%g)^2 )\n = (%g)*( (%g) + (%g) )\n = (%g)\n",
	 exp(-2*x),PI,K[i+1],Gcalc[i+1],K[i-1],Gcalc[i-1],Dl,Gcalc[i+1],Gcalc[i],Gcalc[i-1],Dl, exp(-2*x),PI*(K[i+1]*Gcalc[i+1]-K[i-1]*Gcalc[i-1])/Dl,(Gcalc[i+1]-2*Gcalc[i]+Gcalc[i-1])/(Dl*Dl), dGdt[0] );
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

vpdens = 2PI int_a0^amax G*a*da
vpintegral += G*a*da
vpdens = 2PI*vpintegral

vpdens = 2PI*a0^2 int_0^lmax G*exp(2l)*dl
vpintegral += G*exp(2l)*dl
vpdens = 2PI*a0^2*vpintegral

*/
