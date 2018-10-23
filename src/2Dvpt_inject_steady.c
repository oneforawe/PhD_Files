//  File Comments  //
//=========================================================================//

/* FILENAME: 2Dvpt_inject_steady.c */
/* VERSION: 2 (2011 May 23 - ...)
            Eqns in terms of Gamma (not y), Gamma should no longer fall below its new equilibrium curve */
/* AUTHOR: Andrew Forrester <aforrester@ucla.edu> */
/* DESCRIPTION:
   * This C program calculates time-dependent properties of a 2D system of superfluid (a thin layer of liquid 4He) before and after an instantaneous quench (temperature drop) starting at temperature Ti (Ti/Tkt = iTfrac) at time t=0 and going to temperature Tq (Tq/Tkt = qTfrac) at the next time step.  The system starts in equilibrium and the quench puts the system out of equilibrium so the properties evolve toward a new equilibrium.
   * The properties calculated are superfluid ratios (K, K/K0=sigma_s/sigma_0) and vortex-pair probability density G at various length-scales (l).
*/
/* EXT FILES: none */
/* COMPILE NOTES:
   * To compile, type "g++ -lm 2Dvpt_inject_steady.c" without the quotes; then to run, type "./a.out".
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
const double iTfrac   = 0.10;         //  initial temperature fraction Ti/Tkt
const double qTfrac   = 0.95;         //  quench temperature fraction Tq/Tkt
const double lmax     = 5.0;         //  1  //  maximum of the length scale l (potentially approximate)
const int    lsteps   = 1250;         //  1250, 5000, 100000 //  lsteps = lmax/dl = 10/0.0001 = 100,000  (make sure this is divisible by ldatamax)
const int    lpts     = lsteps+1;     //  from l=0 to l=lmax, inclusive (incl. the boundary-condition-enforcing pnt, there are lpts+1 "l" pnts)
const int    ldatamax = 250;          //  50,100 //  max number of data points recorded per time step
const double linj     = 3.0;          //  injection scale, at which vortices are being injected
const double InjConst = 5.0e0;       //  injection constant - the factor multiplying the Heaviside step function, gamma*t0/(2*PI*a02), where gamma is the factor multiplying the dirac delta source
const double tmax     = 10000.0;      //  20600 //  max unitless time (see if(t==...) below)
const double dt0      = 1.0e-5;       //  the time increment (in units of the "diffusion time")
const double a0       = 1.0;          //  a0 in units of a0 (!)
const double a02      = a0*a0;        //  a0 to the second power
const double a04      = a0*a0*a0*a0;  //  a0 to the fourth power
const double K0c      = 0.747853;     //  critical (T=Tkt) bare superfluid ratio (or "coupling constant")

// Global variables //
int i;
double K[lpts+1], G[lpts+1];
double dGdt[1];



//  Function prototypes  //
//=========================================================================//

// rk4((*f)(),*x,y[],h,n);
// SteadyStateRecRel(x,z,dzdx,n);
// fpKRecRel(x,Kcalc,dKdx,n);
// fpGRecRel(x,Dl,Gcalc);
void rk4(void (*f)(double, double*, double*, unsigned int), double *x, double y[], double h, unsigned n);
void SteadyStateRecRel(double x, double z[2], double dzdx[2], unsigned n);
void fpKRecRel(double x, double Kcalc[1], double dKdx[1], unsigned n);
void fpGRecRel(double x, double Dl, double Gcalc[lsteps+1]);



//  Function definitions  //
//=========================================================================//

main(){
	// Main function definitions //
	int64_t n;
	int j,k;
	double dblsteps=lsteps, dl=lmax/dblsteps;
	double x, l[lpts];
	double K0, G0, G0final;
	double z[2], Kcalc[lpts+1], Gcalc[lpts+1];
	double k1[lpts], k2[lpts], k3[lpts];
	double testG1[lpts+1], testG2[lpts+1];
	double t, dt, ttrig;
	double vpintegral, vpdens;
	FILE *outfile1;  // recording functions versus length-scale l ("vsl"), at selected times
	FILE *outfile2;  // recording functions versus time t ("vst"), at largest length scale lmax
	char *filename1, *filename2;

	// Prepare output files, print identification and values //
	asprintf(&filename1, "2Dvpt_inject_steady_vsl_T_%g_lmax%g_dl%g_linj%g_InjConst%g.dat", iTfrac, lmax, dl, linj, InjConst);
//	asprintf(&filename2, "2DKT_FPinject_SteadyState_vst_T_%g_lmax%g_dl%g_linj%g_InjConst%g_sink.dat", iTfrac, lmax, dl, linj, InjConst);
	outfile1 = fopen(filename1,"w");  //  E.g., "2Dvpt_inject_steady.out"
//	outfile2 = fopen(filename2,"w");  //  E.g., "2Dvpt_inject_steady.out"
	fprintf(outfile1,"# Filename: %s\n", filename1);
//	fprintf(outfile2,"# Filename: %s\n", filename2);
	fprintf(outfile1,"# Source: 2Dvpt_inject_steady.c\n");
//	fprintf(outfile2,"# Source: 2Dvpt_inject_steady.c\n");
	fprintf(outfile1,"# Source version: %s\n", "2 (2011 May 23 - ...)");
//	fprintf(outfile2,"# Source version: %s\n", "2 (2011 May 23 - ...)");
	fprintf(outfile1,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g, B=%21.21g\n", PI,PISQ,PICU,B);
//	fprintf(outfile2,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g, B=%21.21g\n", PI,PISQ,PICU,B);
	fprintf(outfile1,"# Parameter values: iTfrac=%g, qTfrac=%g, lmax=%g, lsteps=%i, (dl=%g), lpts=%i, ldatamax=%i, linj=%g, InjConst=%g, tmax=%g, dt0=%g, a0=%g, a02=%g, a04=%g, K0c=%g\n", iTfrac,qTfrac,lmax,lsteps,dl,lpts,ldatamax,linj,InjConst,tmax,dt0,a0,a02,a04,K0c);
//	fprintf(outfile2,"# Parameter values: iTfrac=%g, qTfrac=%g, lmax=%g, lsteps=%i, (dl=%g), lpts=%i, ldatamax=%i, linj=%g, InjConst=%g, tmax=%g, dt0=%g, a0=%g, a02=%g, a04=%g, K0c=%g\n", iTfrac,qTfrac,lmax,lsteps,dl,lpts,ldatamax,linj,InjConst,tmax,dt0,a0,a02,a04,K0c);

	// Boundary condition enforcement //
	K[lpts] = G[lpts] = 0.0;  // K=G=0 at l=lmax+dl, terminates recursion relations at l=lmax


	///////////////////////////////////////////////
	//                                           //
	//  (n=0) Initial Steady State Calculations  //
	//                                           //
	///////////////////////////////////////////////

	// Initialize time quantities //
	n = 0;
	t = 0.0;
	dt = dt0;  // We won't step in time until the quench.

	// Initialize quantities at smallest length-scale (l=0) //
	x = l[0] = 0.0;
	z[0] = K[0] = K0 = K0c/iTfrac;
	z[1] = G[0] = G0 = exp(-4.0*l[0]-PISQ*K0)/a04;
	vpintegral = 0.0;
	vpintegral += G[0]*exp(2.0*l[0])*dl;

	// Print initial data to screen and output files //
	printf(          "\n# time step n = %ld\tt = %e\t(dt = %e)\n", n,t,dt);
	fprintf(outfile1,"\n# time step n = %ld\tt = %e\t(dt = %e)\n", n,t,dt);
	printf(          "%s\t%s\t%s\t%s\t%s\t%s\n",   "T/Tc","1-T/Tc","l","K","K/K0","G");
	fprintf(outfile1,"# %s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","G");
//	fprintf(outfile2,"# %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "time-step", "t", "T/Tc","1-T/Tc","l","K","K/K0","G","vpdens");
	printf(          "%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[0], K[0], K[0]/K0, G[0]);
	fprintf(outfile1,"%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[0], K[0], K[0]/K0, G[0]);

	// Calculate state-variables at progressively larger length scales, up to coherence length, via Runge-Kutta method (rk4) //
	for(i=1; i<lpts; i++){
		// Step out in length scale l, calculating next l, K and G //
		rk4(SteadyStateRecRel, &x, z, dl, 2);  // equil: rk4 2D (K and G), SteadyStateRecRel
		l[i] = l[i-1] + dl;
		K[i] = z[0];
		G[i] = z[1];
		vpintegral += G[i]*exp(2.0*l[i])*dl;
		// Print at most [ldatamax] (e.g., 50) data points to screen and output file 1 //
		if(i%(lsteps/ldatamax)==0){
			printf(          "%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, G[i]);
			fprintf(outfile1,"%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, G[i]);
		}
	}
	vpdens = 2.0*PI*a02*vpintegral;
	fflush(outfile1);

	// Print out the data at largest length scale to file 2 //
//	fprintf(outfile2,"%ld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", n, t, iTfrac, 1-iTfrac, l[lpts-1], K[lpts-1], K[lpts-1]/K0, G[lpts-1], vpdens);
//	fflush(outfile2);

/*
	//////////////////////////////////////////////////////////////
	//        (After Removing/Changing Injection Rate and/or)   //
	//  (n>=1) After Instantaneous "Quench" (Temperature Drop)  //
	//                                                          //
	//////////////////////////////////////////////////////////////

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
	while(t<=tmax+dt){
		n++;
		t += dt;
		vpintegral = 0.0;
		printf(          "# time step n = %ld\tt = %e\t(dt = %e)\t", n,t,dt);
		//fprintf(outfile1,"# time step n = %ld\tt = %e\t(dt = %e)\t", n,t,dt);

		// Calculate dGdt for all length scales, and update G using third order Runge-Kutta method //
			// nequil: rk3 1D (G), fpGRecRel //
		for(i=1; i<lpts-1; i++){
			fpGRecRel(l[i], dl, G);
			k1[i]     = dGdt[0]*dt;
			testG1[i] = G[i] + 0.5*k1[i];
		}
		k1[0]     = k1[1];                             //*could be questioned
		testG1[0] = G[0] + 0.5*k1[0];                  //*
		k1[lpts-1]     = 2.0*k1[lpts-2] - k1[lpts-3];  //*
		testG1[lpts-1] = G[lpts-1] + 0.5*k1[lpts-1];   //*

		for(i=1; i<lpts-1; i++){
			fpGRecRel(l[i], dl, testG1);
			k2[i]     = dGdt[0]*dt;
			testG2[i] = G[i] - k1[i] + k2[i];
		}
		k2[0]     = k2[1];                                     //*
		testG2[0] = G[0] - k1[0] + k2[0];                      //*
		k2[lpts-1]     = 2.0*k2[lpts-2] - k2[lpts-3];          //*
		testG2[lpts-1] = G[lpts-1] - k1[lpts-1] + k2[lpts-1];  //*

		for(i=1; i<lpts-1; i++){
			fpGRecRel(l[i], dl, testG2);
			k3[i] = dGdt[0]*dt;
			G[i]  = G[i] + k1[i]/6.0 + 2.0*k2[i]/3.0 + k3[i]/6.0;
			vpintegral += G[i]*exp(2*l[i])*dl;
		}
		// G[0] = G0final; // might want to use this bndy condition
		G[0] = 2.0*G[1] - G[2];                       //*
		if(G[0]<=G0final)  G[0] = G0final;            //*
		G[lpts-1] = G[lpts-2]/(1+2*PI*K[lpts-1]*dl);  //* WHY?
		// The line above is based on the (false?) equilibrium relation dG/dl = -2piKG.
		// Could instead use: G[lpts-1] = G[lpts-2] - 2*PI*K[lpts-2]*dl);
		//                or: G[lpts-1] = G[lpts-2] + G[lpts-2]*(4.0-4.0*l[lpts-2]-2*PI*K[lpts-2])*dl;
		//                or: G[lpts-1] = G[lpts-2]/(1.0-(4.0-4.0*l[lpts-1]-2*PI*K[lpts-1])*dl);
		// Shouldn't the above be based on... something else (non-equil)?
		vpintegral += G[0]*exp(2*l[0])*dl;
		vpintegral += G[lpts-1]*exp(2*l[lpts-1])*dl;
		vpdens = 2*PI*a02*vpintegral;

		// Calculate K at progressively larger length scales, up to coherence length, via Runge-Kutta method (rk4 K) //
		x = 0.0;
		Kcalc[0] = K[0];
		for(i=1; i<lpts; i++){
			rk4(fpKRecRel, &x, Kcalc, dl, 1);  // nequil: rk4 1D (K), fpKRecRel
			K[i] = Kcalc[0];
		}

		// At certain times, print out the results to file 1 //
		if( (t>=0.01&&j==1) || (t>=0.1&&j==2) || (t>=1&&j==3) || (t>=10&&j==4) || (t>=100&&j==5) || (t>=1000&&j==6) || (t>=10000&&j==7) ){
			fprintf(outfile1,"\n# time step n = %ld\tt = %e\t(dt = %e)\n", n,t,dt);
			fprintf(outfile1,"# %s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","G");
			for(i=0; i<lpts; i++){
				// Print at most [ldatamax] (e.g., 1000) data points to screen and output file //
				if(i%(lsteps/ldatamax)==0){
					printf(          "%g\t%g\t%g\t%g\t%g\t%g\n", qTfrac, 1-qTfrac, l[i], K[i], K[i]/K0, G[i]);
					fprintf(outfile1,"%g\t%g\t%g\t%g\t%g\t%g\n", qTfrac, 1-qTfrac, l[i], K[i], K[i]/K0, G[i]);
				}
			}
			fflush(outfile1);
			j++;
		}

		// Print out the data at the largest length-scale to file 2, selecting data at exponentially sparser time intervals //
		if(t>ttrig){
			fprintf(outfile2,"%ld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", n, t, qTfrac, 1-qTfrac, l[lpts-1], K[lpts-1], K[lpts-1]/K0, G[lpts-1], vpdens);
			fflush(outfile2);
			ttrig += dt0*pow(10.0,0.05*k);
			//ttrig += dt0*exp(0.05*k);
			k++;
		}

		// Bug detection //
		printf(          "%g\t%g\t%g\n", l[lpts-1], K[lpts-1]/K0, G[lpts-1]);
		//fprintf(outfile1,"%g\t%g\t%g\n", l[lpts-1], KG[lpts-1].c[0]/K0, KG[lpts-1].c[1]);
	}
*/
	fclose(outfile1);
//	fclose(outfile2);
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



// equilibrium (K,G) recursion relations (No Dirac-Delta Sink and all outflow at a0) //
void SteadyStateRecRel(double x, double z[2], double dzdx[2], unsigned n){
	dzdx[0] = -B*a04*exp(4.0*x)*z[1]*z[0]*z[0];
	if(x<linj)
		dzdx[1] = -2.0*PI*z[0]*z[1] + InjConst;
	else
		dzdx[1] = -2.0*PI*z[0]*z[1];
}



/*
// equivalent results using this algorithm:
// equilibrium (K,G) recursion relations (w/ Dirac-Delta Sink and zero-flow at a0) //
void SteadyStateRecRel(double x, double z[2], double dzdx[2], unsigned n){
	dzdx[0] = -B*a04*exp(4.0*x)*z[1]*z[0]*z[0];
	if(x==0||x>linj)
		dzdx[1] = -2.0*PI*z[0]*z[1];
	else
		dzdx[1] = -2.0*PI*z[0]*z[1] + InjConst;
}
*/



// Fokker-Planck K recursion relations //
void fpKRecRel(double x, double Kcalc[1], double dKdx[1], unsigned n){
	dKdx[0] = -B*a04*exp(4*x)*G[i-1]*Kcalc[0]*Kcalc[0];
}



// Fokker-Planck G recursion relations //
void fpGRecRel(double x, double Dl, double Gcalc[lsteps+1]){
	// My version of the equation:
	dGdt[0] = exp(-2.0*x)*( PI*(K[i+1]*Gcalc[i+1]-K[i-1]*Gcalc[i-1])/Dl + (Gcalc[i+1]-2.0*Gcalc[i]+Gcalc[i-1])/(Dl*Dl) );
	// Han-Ching Chu's version, which is mine divided by 2*PI (AND k is multiplied with the unmodified g, as opposed to gin):
	//dGdt[0] = exp(-2.0*x)*( (K[i+1]*Gcalc[i+1]-K[i-1]*Gcalc[i-1])/2.0/Dl + (1/2.0/PI)*(Gcalc[i+1]-2.0*Gcalc[i]+Gcalc[i-1])/(Dl*Dl) );
	//gstep[0]= exp(-2.0*x)*( (k[j+1]*g[j+1]-k[j-1]*g[j-1])/2.0/stepsize + (1/2.0/PI)*(gin[j+1]-2.0*gin[j]+gin[j-1])/(stepsize*stepsize) );
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
