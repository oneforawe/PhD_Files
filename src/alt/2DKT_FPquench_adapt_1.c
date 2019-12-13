//  File Comments  //
//=========================================================================//

/* FILENAME: 2DKT_FPquench_adapt_1.c */
/* VERSION: 1 (2011 Jan 02)
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
const double lmax      = 10.0;      //  1  //  maximum of the length scale l (potentially approximate)
const int    lsteps    = 1250;      //  5000 // 100000 //  lsteps = lmax/dl = 10/0.0001 = 100,000  (make sure this is divisible by ldatamax)
const int    lpts      = lsteps+1;  //  from l=0 to l=lmax, inclusive (incl. the boundary-condition-enforcing pnt, there are lpts+1 "l" pnts)
const int    ldatamax  = 50;        //  maximum number of l data points recorded per time step (beyond the l=0 data point)
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
const double PSHRINK   = -0.25;
const double PGROW     = -0.20;
const double MAXFACTOR = 3.0;
const double SAFETY    = 0.9;

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

// Global variables //
int i;
double l[lpts], dl[lpts], K[lpts+1], G[lpts+1];
double dGdt[1];



//  Function prototypes  //
//=========================================================================//
// my_rkf45((*derivs)(),*x,h,y[],n,dre,*htemp,ytemp[],yerr[]);
// EqRecRel(x,z[2],dzdx[2],n);
// fpKRecRel(x,Kcalc[1],dKdx[1],n);
// fpGRecRel(x,dx,Gcalc[lsteps+1]);
void my_rkf45( void (*derivs)(double, double*, double*, unsigned),
	       double *x, double h, double y[], unsigned n, double dre, double *htemp, double ytemp[], double yerr[] );
void EqRecRel(double x, double z[2], double dzdx[2], unsigned n);
void fpKRecRel(double x, double Kcalc[1], double dKdx[1], unsigned n);
void fpGRecRel(double x, double dx, double Gcalc[lsteps+1]);



//  Function definitions  //
//=========================================================================//

main(){
	// Main function definitions //
	int64_t n;
	int imax,accuracy=0,j,k=1;
	double t, dt;
	double dblmax=lmax, dblsteps=lsteps, dbldatamax=ldatamax;
	double dl0=dblmax/dblsteps, dataDl=lmax/dbldatamax;
	double dltemp;
	double K0, y0;
	double z[2], ztemp[2], zerr[2];
	double Kcalc[1], Kcalctemp[1], Kcalcerr[1];
	double Ktemp[lpts], Kerr[lpts], Gerr[lpts];
	double KerrTot, GerrTot;
	double errfrac;
	double k1[lpts], k2[lpts], k3[lpts];
	double testG1[lpts+1], testG2[lpts+1];
	FILE *outfile;
	char filename[100];

	// Prepare output file, print identification and values //
	sprintf(filename, "2DKT_FPquench_adapt_lmax%g_dl0%g_tmax%g_dt0%g_DRE%g.dat", lmax, dl0, tmax, dt0, DRE);
	outfile = fopen(filename,"w");  //  E.g., "2DKT_FPquench_adapt.out"
	fprintf(outfile,"# Filename: %s\n", filename);
	fprintf(outfile,"# Source: 2DKT_FPquench_adapt.c\n");
	fprintf(outfile,"# Source version: %s\n", "1 (2010 12 09)");
	fprintf(outfile,"# Parameter values: lmax=%g, lsteps=%i, (dl0=%g), lpts=%i, ldatamax=%i, iTfrac=%g, qTfrac=%g, tmax=%g, dt0=%g, a0=%g, a04=%g, K0c=%g, DRE=%g, PSHRINK=%g, PGROW=%g, MAXFACTOR=%g\n", lmax,lsteps,dl0,lpts,ldatamax,iTfrac,qTfrac,tmax,dt0,a0,a04,K0c,DRE,PSHRINK,PGROW,MAXFACTOR);

	// Boundary condition enforcement //
////	K[lpts] = G[lpts] = 0;  // K=G=0 at l=lmax+dl, terminates recursion relations at l=lmax

	//////////////////////////////////////////////
	//                                          //
	//  (n=0) Initial Equilibrium Calculations  //
	//                                          //
	//////////////////////////////////////////////

	// Initialize time quantities //
	n = 0;
	t = 0.0;
	dt = dt0;

	// Initialize temperature-dependent quantities at smallest length-scale (l=0) //
	l[0] = 0.0;
	dl[0] = dl0;
	K0 = K0c/iTfrac;
	y0 = exp(-PISQ*K0/2);
	z[0] = K0;
	z[1] = y0;
	K[0] = z[0];
	G[0] = z[1]*z[1]*exp(-4.0*l[0])/a04;  // G (not = 0.0)

	// Print headings and initial data to screen and output file //
	fprintf(outfile,"\n# time step n = %ld\tt = %e\t(dt = %e)\n", n,t,dt);
	printf(         "%s\t%s\t%s\t%s\t%s\t%s\n",   "T/Tc","1-T/Tc","l","K","K/K0","G");
	fprintf(outfile,"# %s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","G");
	printf(         "%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[0], K[0], K[0]/K0, G[0]);
	fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[0], K[0], K[0]/K0, G[0]);

	// Calculate state-variables at progressively larger length scales, up to coherence length, via Runge-Kutta method //
	i = 0; // for stepping l
	j = 1; // for selecting data to print to file
	while(l[i]<lmax){
		i++;
		imax = i;
		//printf("i = %i\n",i);
		// Step out in length scale l, calculating next K, y, l, and G //
		accuracy = 0;
		while(accuracy==0){
			my_rkf45(EqRecRel, &l[i-1], dl[i-1], z, 2, DRE, &dltemp, ztemp, zerr);
			// equil: rkf6 2D (K and y), EqRecRel, adaptive step to desired relative error DRE
			//printf("%s\n","made it out of rk!");
			//printf("\tdltemp = %g\n",dltemp);
			l[i] = l[i-1]+dltemp;
			G[i] = ztemp[1]*ztemp[1]*exp(-4.0*l[i])/a04;
			//printf("\tG[%i] = %g\n",i,G[i]);
			Gerr[i] = zerr[1]*2*ztemp[1]*exp(-4.0*l[i])/a04;
			//printf("\tGerr[%i] = %g\n",i,Gerr[i]);
			// check error (in G), redo (w/ smaller increment) or proceed (w/o changing increment) //
			errfrac = fabs( Gerr[i]/(DRE*G[i]) );  // (calculated error)/(desired error)
			//printf("\terrfrac = %g\n",errfrac);
			if(errfrac>1.0){  // errfrac<1.0 means calculated error is too big
				accuracy = 0;
				dl[i-1] = SAFETY*dltemp*pow(errfrac,0.2);
				// Failure print out //
				//printf("G accuracy Fail: %g\t%g\t%g\t%g\t%g\n", l[i-1], dltemp, l[i], ztemp[0], G[i]);
			}  // otherwise, use the dl already defined by my_rkf6
			else{
				accuracy = 1;
				dl[i-1] = dltemp;
				dl[i] = dltemp;
				z[0] = ztemp[0];
				K[i] = z[0];
				Kerr[i] = zerr[0];
				z[1] = ztemp[1];
				//printf("G accuracy Success: %g\t%g\t%g\t%g\t%g\n", l[i-1], dltemp, l[i], ztemp[0], G[i]);
			}
		}
		// Print at most [ldatamax] (e.g., 50) data points to screen and output file //
		if(l[i]>=j*dataDl){
			while(l[i]>=j*dataDl)  j++;
			printf(         "%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, G[i]);
			fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, G[i]);
		}
	}

	// Calculate and print total error //
	for(i=0;i<=imax;i++){
		GerrTot += Gerr[i];
		KerrTot += Kerr[i];
	}
	fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[imax], (K[imax]+KerrTot), (K[i]+KerrTot)/K0, (G[i]+GerrTot));
	fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[imax], (K[imax]-KerrTot), (K[i]-KerrTot)/K0, (G[i]-GerrTot));
	fprintf(outfile,"# Total error in K = %g\n", KerrTot);
	fprintf(outfile,"# Total error in G = %g\n", GerrTot);
	fprintf(outfile,"# By the way, imax = %i\n", imax);
	fflush(outfile);
/*
	//////////////////////////////////////////////////////////////
	//                                                          //
	//  (n>=1) After Instantaneous "Quench" (Temperature Drop)  //
	//                                                          //
	//////////////////////////////////////////////////////////////

	// Initialize temperature-dependent quantities at smallest length-scale (l=0) //
	l[0] = 0.0;
	dl[0] = dl0;
	K0 = K0c/qTfrac;
	//y0 = exp(-PISQ*K0/2);
	Kcalc[0] = K0;
	K[0] = Kcalc[0];
	//G[0] = y0*y0*exp(-4.0*l[0])/a04;  // G (not = 0.0)


	// Calculate K at progressively larger length scales, up to coherence length, via Runge-Kutta method (derkK) //
	while(l<lmax){
		my_rkf45(fpKRecRel, &l, Kcalc, dl, 1);  // nequil: rkf6 1D (K), fpKRecRel
		K = Kcalc[0];
		// check error (in K), redo (w/ smaller increment) or proceed (setting larger increment) //
	}

	i = 0; // for stepping l
	j = 1; // for selecting data to print to file
	while(l[i]<lmax){
		i++;
		// Step out in length scale l, calculating next K, and l //
		my_rkf45(fpKRecRel, &l[i-1], dl[i-1], Kcalc, 1, DRE, &dltemp, Kcalctemp, Kcalcerr);
		// equil: rkf6 2D (K and y), EqRecRel, adaptive step to desired relative error DRE
		dl[i-1] = dltemp;
		l[i] = l[i-1]+dl[i-1];
		dl[i] = dltemp;
		Kcalc[0] = Kcalctemp[0];
		K[i] = Kcalc[0];
		Kerr[i] = zerr[0];
		z[1] = ztemp[1];
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
			//derk(fpKRecRel, &x, Kcalc, dl, 1);  // nequil: rk4 1D (K), fpKRecRel
			my_rkf6(fpKRecRel, &x, Kcalc, dl, 1);  // nequil: rkf6 1D (K), fpKRecRel
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
*/
	fclose(outfile);
}



void my_rkf45(void (*derivs)(double, double*, double*, unsigned), double *x, double h, double y[], unsigned n, double dre, double *htemp, double ytemp[], double yerr[]){
	int i,accurate=0,check=0;
	double *k[6],*ycalc,maxerr,temperr,factor;

	for(i=0;i<6;i++)  k[i] = (double *)calloc(n,sizeof(double));  // the k's are slopes, dy/dx
	ycalc = (double *)calloc(n,sizeof(double));

	*htemp = h;
	while(accurate==0){
//	for(check=0;check<20;check++){
		maxerr=1e-300;
		derivs(*x,y,k[0],n);           // get k1=k[0]
		for(i=0;i<n;i++)  ycalc[i] = y[i]+*htemp*(b21*k[0][i]);
		derivs(*x+*htemp*a2,ycalc,k[1],n);  // get k2=k[1]
		for(i=0;i<n;i++)  ycalc[i] = y[i]+*htemp*(b31*k[0][i]+b32*k[1][i]);
		derivs(*x+*htemp*a3,ycalc,k[2],n);  // get k3=k[2]
		for(i=0;i<n;i++)  ycalc[i] = y[i]+*htemp*(b41*k[0][i]+b42*k[1][i]+b43*k[2][i]);
		derivs(*x+*htemp*a4,ycalc,k[3],n);  // get k4=k[3]
		for(i=0;i<n;i++)  ycalc[i] = y[i]+*htemp*(b51*k[0][i]+b52*k[1][i]+b53*k[2][i]+b54*k[3][i]);
		derivs(*x+*htemp*a5,ycalc,k[4],n);  // get k5=k[4]
		for(i=0;i<n;i++)  ycalc[i] = y[i]+*htemp*(b61*k[0][i]+b62*k[1][i]+b63*k[2][i]+b64*k[3][i]+b65*k[4][i]);
		derivs(*x+*htemp*a6,ycalc,k[5],n);  // get k6=k[5]

		for(i=0;i<n;i++){
			ytemp[i] = y[i] + *htemp*(c1*k[0][i]+c2*k[1][i]+c3*k[2][i]+c4*k[3][i]+c5*k[4][i]+c6*k[5][i]);
			//printf("\tytemp[%i] = %g\t",i,ytemp[i]);
			yerr[i] = *htemp*(dc1*k[0][i]+dc2*k[1][i]+dc3*k[2][i]+dc4*k[3][i]+dc5*k[4][i]+dc6*k[5][i]);
			//printf("yerr[%i] = %g\t",i,yerr[i]);
			temperr = fabs( yerr[i]/(dre*ytemp[i]) );  // (calculated error)/(desired error)
			//printf("temperr = %g\n",temperr);
			if(temperr>maxerr)  maxerr = temperr;
		}
		//printf("\tmaxerr = %g\n",maxerr);
		if(maxerr>1.0){  // calculated error is too big
			accurate = 0;
			// Failure print out //
			//printf("rkfail: %g\t%g\t%g\t%g\t%g\t%g\t%g\t%i\n", *x, *htemp, ytemp[0], yerr[0], ytemp[1], yerr[1], maxerr, accurate);
			*htemp *= SAFETY*pow(maxerr,PSHRINK);
		}
		else{  // calculated error small enough (or too small, so grow h!)
			accurate = 1;
			factor = SAFETY*pow(maxerr,PGROW);
			if(factor<MAXFACTOR) *htemp *= factor;
			else *htemp *= MAXFACTOR;
		}
		//else  for(i=0;i<n;i++)  y[i] = ytemp[i];
		//htemp *= pow(minfrac,0.2);
		// what's going on? //
		//printf("rk: accurate = %i\thtemp = %g\n", accurate, *htemp);
	}

	for(i=0;i<6;i++)  free((char *)k[i]);
	free((char *)ycalc);
}



void EqRecRel(double x, double z[2], double dzdx[2], unsigned n){
	dzdx[0] = -B*z[1]*z[1]*z[0]*z[0];
	dzdx[1] = (2-PI*z[0])*z[1];
}



void fpKRecRel(double x, double Kcalc[1], double dKdx[1], unsigned n){
	dKdx[0] = -4*PI*PI*PI*G[i-1]*Kcalc[0]*Kcalc[0]*exp(4.0*x);
}



void fpGRecRel(double x, double dx, double Gcalc[lsteps+1]){
	dGdt[0] = exp(-2.0*x)*((1/2.0/PI)*(Gcalc[i+1] - 2.0*Gcalc[i] + Gcalc[i-1])/(dx*dx) + (K[i+1]*Gcalc[i+1] - K[i-1]*Gcalc[i-1])/2.0/dx);
	//dzdt = exp(-2*l)*( (zhi[0]*zhi[1]-zlo[0]*zlo[1])/2.0/Dl + (1/2.0/PI)*(zhi[1]-2*zmd[1]+zlo[1])/(Dl*Dl) );
	//dzdt = exp(-2*l)*( PI*(zhi[0]*zhi[1]-zlo[0]*zlo[1])/Dl + (zhi[1]-2*zmd[1]+zlo[1])/(Dl*Dl) );
}



/*  Program Notes  */
/*=========================================================================*/
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


