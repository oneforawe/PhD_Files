/*  File Comments  */
/*=========================================================================*/

// For some reason, this version has a weird glitchy drop in the Gamma versus a/a0 data from the first data point to the second data point.

/* FILENAME: 2DKT_FPquench.c */
/* VERSION: 1 (2010 12 09)
/* AUTHOR: Andrew Forrester <aforrester@ucla.edu> */
/* DESCRIPTION:
   * This C program calculates time-dependent properties of a 2D system of superfluid (a thin layer of liquid 4He) before and after an instantaneous quench (temperature drop) starting at temperature Ti (Ti/Tkt = iTfrac) at time t=0 and going to temperature Tq (Tq/Tkt = qTfrac) at the next time step.  The system starts in equilibrium and the quench puts the system out of equilibrium so the properties evolve toward a new equilibrium.
   * The properties calculated are superfluid ratios (K, K/K0 = rho/rho0) and vortex-pair probability density G at various length-scales (l).
*/
/* EXT FILES: none */
/* COMPILE NOTES:
   * To compile, type "g++ -lm 2DKT_FPquench.c" without the quotes; then to run, type "./a.out".
   * NOTE: Be sure...
*/

// WORK ON 2DKT_FPquench_adapt.c!!!

/*  Function Preparation  */
/*=========================================================================*/

/* Standard routine header files */
#include <stdio.h>
#include <math.h>
#include <stdint.h>

/* Constants */
const double PI   = 3.14159265358979323846;
const double PISQ = 9.86960440108935799230;
const double PICU = 31.00627668029981620634;
const double B    = 4.0*PICU;
const int    N    = 2;

/* Program parameters/inputs and data type definitions */
const int    lmax     = 10;        //  1  //  maximum of the length scale l (potentially approximate)
const int    lsteps   = 1250;      //  5000 // 100000 //  lsteps = lmax/dl = 10/0.0001 = 100,000  (make sure this is divisible by ldatamax)
const int    lpts     = lsteps+1;  //  from l=0 to l=lmax, inclusive (incl. the boundary-condition-enforcing pnt, there are lpts+1 "l" pnts)
const int    ldatamax = 50;        //  100 //  max number of data points recorded per time step
const double iTfrac   = 0.95;      //  initial temperature fraction Ti/Tkt
const double qTfrac   = 0.10;      //  quench temperature fraction Tq/Tkt
const double tmax     = 10000.0;   //  20600 //  max unitless time (see if(t==...) below)
const double dt0      = 1.8e-4;    //  the time increment (in units of the "diffusion time")
const double a0       = 1.0;       //  a0 in units of a0 (!)
const double a04      = 1.0;       //  a0 to the fourth power
const double K0c      = 0.747853;  //  critical (T=Tkt) bare superfluid ratio (or "coupling constant")

typedef struct {
	double c[N];  //  c = "component/choose"
} RETARRAY;

double l[lpts], dGdt[lpts];
RETARRAY Ky[lpts], KG[lpts+1];

double testDGDt[lpts+1][4], testdt;
double k1[lpts+1], k2[lpts+1], k3[lpts+1], k4[lpts+1];
RETARRAY testKG[lpts+1][3];



/*  Function prototypes  */
/*=========================================================================*/
//  derk(func(),CalcKy,Calcl,Dl,n)
//  derkK(func(),CalcKG,Calcl,Dl,n)
//  EqRecRel(l,z)
//  fpKRecRel(l,z)
//  fpGRecRel(l,Dl,zlo,zmd,zhi);
RETARRAY derk(RETARRAY (*func)(double l, double z[]), RETARRAY CalcKy, double Calcl, double Dl, unsigned n);
double   derkK(double (*func)(double l, double z[]), RETARRAY CalcKG, double Calcl, double Dl);
RETARRAY EqRecRel(double l, double z[]);
double   fpKRecRel(double l, double z[]);
double   fpGRecRel(double l, double Dl, double zlo[], double zmd[], double zhi[]);



/*  Function definitions  */
/*=========================================================================*/

main(){
	/* Main function definitions */
	int64_t n;
	int i,j,k=1;
	double dblmax=lmax, dblsteps=lsteps;
	double K0, y0, dl=dblmax/dblsteps;
	double t, dt, dtLimit;
	FILE *outfile;
	char filename[100];

	/* Prepare output file, print identification and values */
	sprintf(filename, "2DKT_FPquench_lmax%i_dl%g_tmax%g_dt0%g.dat", lmax, dl, tmax, dt0);
	outfile = fopen(filename,"w");  //  E.g., "2DKT_FPquench.out"
	fprintf(outfile,"# Filename: %s\n", filename);
	fprintf(outfile,"# Source: 2DKT_FPquench.c\n");
	fprintf(outfile,"# Source version: %s\n", "1 (2010 12 09)");
	fprintf(outfile,"# Parameter values: lmax=%i, lsteps=%i, (dl=%g), lpts=%i, ldatamax=%i, iTfrac=%g, qTfrac=%g, tmax=%g, dt0=%g, a0=%g, a04=%g, K0c=%g\n", lmax,lsteps,dl,lpts,ldatamax,iTfrac,qTfrac,tmax,dt0,a0,a04,K0c);

	/* Boundary condition enforcement */
	KG[lpts].c[0] = KG[lpts].c[1] = 0;  // K=G=0 at l=lmax+dl, terminates recursion relations at l=lmax

	/* Define the timestep size */
	// dtLimit = 2.5*(dl*dl)/sqrt((4.0/PI/PI) + K0*dl*K0*dl);
	dt = dt0;  //0.9*dtLimit;

	/********************************************/
	/*                                          */
	/*  (n=0) Initial Equilibrium Calculations  */
	/*                                          */
	/********************************************/

	/* Initialize */
	n = 0;
	t = 0.0;
	//dt = 0.0;  //  no time steps yet
	l[0] = 0.0;
	K0 = K0c/iTfrac;
	y0 = exp(-PISQ*K0/2);
	Ky[0].c[0] = K0;  //  K
	KG[0].c[0] = K0;  //  K
	Ky[0].c[1] = y0;  //  y
	KG[0].c[1] = exp(-4*l[0])*y0*y0/a04;  // G (not = 0.0)

	/* Print initial data to screen and output file */
	fprintf(outfile,"\n# time step n = %ld\tt = %e\t(dt = %e)\n", n,t,dt);
	printf(         "%s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","G");
	fprintf(outfile,"# %s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","G");
	printf(         "%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[0], KG[0].c[0], KG[0].c[0]/K0, KG[0].c[1]);
	fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[0], KG[0].c[0], KG[0].c[0]/K0, KG[0].c[1]);

	/* Calculate state-variables at progressively larger length scales, up to coherence length, via Runge-Kutta method (derk) */
	for(i=1; i<lpts; i++){
		/* Calculate next Ky and KG state, incrementing length scale */
		Ky[i] = derk(EqRecRel, Ky[i-1], l[i-1], dl, N);
		l[i] = l[i-1] + dl;
		KG[i].c[0] = Ky[i].c[0];
		KG[i].c[1] = exp(-4*l[i])*Ky[i].c[1]*Ky[i].c[1]/a04;

		/* Print at most [ldatamax] (e.g., 50) data points to screen and output file */
		if(i%(lsteps/ldatamax)==0){
			printf(         "%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[i], KG[i].c[0], KG[i].c[0]/K0, KG[i].c[1]);
			fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[i], KG[i].c[0], KG[i].c[0]/K0, KG[i].c[1]);
		}
	}
	fflush(outfile);

	/************************************************************/
	/*                                                          */
	/*  (n>=1) After Instantaneous "Quench" (Temperature Drop)  */
	/*                                                          */
	/************************************************************/

	/* Initialize */
	// l[0] = 0.0;  //  (already done)
	K0 = K0c/qTfrac;
	//y0 = exp(-PISQ*K0/2.0);
	KG[0].c[0] = K0;  //  K
	//KG[0].c[1] = exp(-4.0*l[0])*y0*y0/a04;  // G (not = 0.0)  //  (correct location?)

	/* Calculate K at progressively larger length scales, up to coherence length, via Runge-Kutta method (derkK) */
	for(i=1; i<lpts; i++)
		KG[i].c[0] = derkK(fpKRecRel, KG[i-1], l[i-1], dl);

	/* Advance in time */
	while(t<=tmax+dt){
		n++;
		t += dt;
		printf(         "# time step n = %ld\tt = %e\t(dt = %e)\t", n,t,dt);
		//fprintf(outfile,"# time step n = %ld\tt = %e\t(dt = %e)\t", n,t,dt);

		/* Calculate dGdt for all length scales, and update G */
/*		dGdt = derkG(fpGRecRel, KG, l, dl, dt);       */
/*		for(j=0; j<4; j++){
			testDGDt[0][j]    = 0;
			testDGDt[lpts][j] = 0;
		}
		for(i=1; i<lpts-1; i++)
			testDGDt[i][0] = fpGRecRel(l[i], dl, KG[i-1].c, KG[i].c, KG[i+1].c);
		for(j=0; j<3; j++){
			testdt = dt*((j+2)/2)/2;  // This algorithm, using '((...)/2)/2', creates a desired rounding "error" that produces the sequence 0.0, 0.5, 0.5, 1.0, instead of 0.0, 0.5, 0.75, 1.0.
			for(i=0; i<lpts; i++){
				testKG[i][j].c[0] = KG[i].c[0];
				testKG[i][j].c[1] = KG[i].c[1] + testDGDt[i][j]*testdt;
			}
			for(i=1; i<lpts-1; i++)
				testDGDt[i][j+1] = fpGRecRel(l[i], dl, testKG[i-1][j].c, testKG[i][j].c, testKG[i+1][j].c);
		}
		for(i=1; i<lpts; i++)
			dGdt[i] = (testDGDt[i][0]+2*(testDGDt[i][1]+testDGDt[i][2])+testDGDt[i][3])/6;  */

		/* Step G in time, for all length scales
		for(i=1; i<lpts; i++)
			KG[i].c[1] += dGdt[i]*dt; */


/*		for(j=0; j<4; j++){
			testDGDt[0][j]    = 0;
			testDGDt[lpts][j] = 0;
		}
 
		for(i=1; i<lpts-1; i++)
			testDGDt[i][0]    = fpGRecRel(l[i], dl, KG[i-1].c, KG[i].c, KG[i+1].c);
		for(i=0; i<lpts; i++){
			testKG[i][0].c[0] = KG[i].c[0];
			k1[i]             = testDGDt[i][0]*dt;
			testKG[i][0].c[1] = KG[i].c[1] + k1[i]/2.0;
		}

		for(i=1; i<lpts-1; i++)
			testDGDt[i][1]    = fpGRecRel(l[i], dl, testKG[i-1][0].c, testKG[i][0].c, testKG[i+1][0].c);
		for(i=0; i<lpts; i++){
			testKG[i][1].c[0] = KG[i].c[0];
			k2[i]             = testDGDt[i][1]*dt;
			testKG[i][1].c[1] = KG[i].c[1] - k1[i] + k2[i];
		}

		for(i=1; i<lpts-1; i++)
			testDGDt[i][2]    = fpGRecRel(l[i], dl, testKG[i-1][1].c, testKG[i][1].c, testKG[i+1][1].c);
		for(i=0; i<lpts; i++){
			testKG[i][2].c[0] = KG[i].c[0];
			k3[i]             = testDGDt[i][2]*dt;
			testKG[i][2].c[1] = KG[i].c[1] + k1[i]/6.0 + 2.0*k2[i]/3.0 + k3[i]/6.0;
		}

		for(i=1; i<lpts-1; i++)
			testDGDt[i][3]    = fpGRecRel(l[i], dl, testKG[i-1][2].c, testKG[i][2].c, testKG[i+1][2].c);
		for(i=0; i<lpts; i++){
			k4[i]             = testDGDt[i][3]*dt;
			KG[i].c[1]        = KG[i].c[1] + k1[i]/6.0 + k2[i]/3.0 + k3[i]/3.0 + k4[i]/6.0;
		}  */


/*  */		for(i=1; i<lpts-1; i++){
			testDGDt[i][0]    = fpGRecRel(l[i], dl, KG[i-1].c, KG[i].c, KG[i+1].c);
			k1[i]             = testDGDt[i][0]*dt;
			testKG[i][0].c[0] = KG[i].c[0];
			testKG[i][0].c[1] = KG[i].c[1] + k1[i]/2.0;
		}
		testKG[0][0].c[1]      = testKG[1][0].c[1];
		testKG[lpts-1][0].c[1] = testKG[lpts-2][0].c[1];

		for(i=1; i<lpts-1; i++){
			testDGDt[i][1]    = fpGRecRel(l[i], dl, testKG[i-1][0].c, testKG[i][0].c, testKG[i+1][0].c);
			k2[i]             = testDGDt[i][1]*dt;
			testKG[i][1].c[0] = KG[i].c[0];
			testKG[i][1].c[1] = KG[i].c[1] - k1[i] + k2[i];
		}
		testKG[0][1].c[1]      = testKG[1][1].c[1];
		testKG[lpts-1][1].c[1] = testKG[lpts-2][1].c[1];

		for(i=1; i<lpts-1; i++){
			testDGDt[i][2]    = fpGRecRel(l[i], dl, testKG[i-1][1].c, testKG[i][1].c, testKG[i+1][1].c);
			k3[i]             = testDGDt[i][2]*dt;
/*			testKG[i][2].c[0] = KG[i].c[0];
			testKG[i][2].c[1] = KG[i].c[1] + k1[i]/6.0 + 2.0*k2[i]/3.0 + k3[i]/6.0;  */
			KG[i].c[1]        = KG[i].c[1] + k1[i]/6.0 + 2.0*k2[i]/3.0 + k3[i]/6.0;  //  update G here and below**
		}
/*		testKG[0][2].c[1]      = testKG[1][2].c[1];
		testKG[lpts-1][2].c[1] = testKG[lpts-2][2].c[1];

		for(i=1; i<lpts-1; i++){
			testDGDt[i][3]    = fpGRecRel(l[i], dl, testKG[i-1][2].c, testKG[i][2].c, testKG[i+1][2].c);
			k4[i]             = testDGDt[i][3]*dt;
			KG[i].c[1]        = KG[i].c[1] + k1[i]/6.0 + k2[i]/3.0 + k3[i]/3.0 + k4[i]/6.0;
		}  */
		KG[0].c[1]      = KG[1].c[1];       //  **here
		KG[lpts-1].c[1] = KG[lpts-2].c[1];  //  **and here


		/* Initialize */
		// l[0] = 0.0;  //  (already done)
		K0 = K0c/qTfrac;
		KG[0].c[0] = K0;  //  K
		//fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\n", qTfrac, 1-qTfrac, l[0], KG[0].c[0], KG[0].c[0]/K0, KG[0].c[1]);
		//printf(         "%g\t%g\t%g\t%g\t%g\t%g\n", qTfrac, 1-qTfrac, l[0], KG[0].c[0], KG[0].c[0]/K0, KG[0].c[1]);

		/* Calculate K at progressively larger length scales, up to coherence length, via Runge-Kutta method (derk1) */
		for(i=1; i<lpts; i++){
			KG[i].c[0] = derkK(fpKRecRel, KG[i-1], l[i-1], dl);
			//if(i%(lsteps/ldatamax)==0){
			//	fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\n", qTfrac, 1-qTfrac, l[i], KG[i].c[0], KG[i].c[0]/K0, KG[i].c[1]);
			//	printf(         "%g\t%g\t%g\t%g\t%g\t%g\n", qTfrac, 1-qTfrac, l[i], KG[i].c[0], KG[i].c[0]/K0, KG[i].c[1]);
			//}
		}

		//fflush(outfile);

		/* At certain times, print out the results */
		if( (t>=0.1&&k==1) || (t>=1&&k==2) || (t>=10&&k==3) || (t>=100&&k==4) || (t>=1000&&k==5) || (t>=10000&&k==6) ){
			fprintf(outfile,"\n# time step n = %ld\tt = %e\t(dt = %e)\n", n,t,dt);
			fprintf(outfile,"# %s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","G");
			for(i=0; i<lpts; i++){
				/* Print at most [ldatamax] (e.g., 1000) data points to screen and output file */
				if(i%(lsteps/ldatamax)==0){
					fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\n", qTfrac, 1-qTfrac, l[i], KG[i].c[0], KG[i].c[0]/K0, KG[i].c[1]);
					printf(         "%g\t%g\t%g\t%g\t%g\t%g\n", qTfrac, 1-qTfrac, l[i], KG[i].c[0], KG[i].c[0]/K0, KG[i].c[1]);
				}
			}
			fflush(outfile);
			k++;
		}

		printf(         "%g\t%g\t%g\n", l[lpts-1], KG[lpts-1].c[0]/K0, KG[lpts-1].c[1]);
		//fprintf(outfile,"%g\t%g\t%g\n", l[lpts-1], KG[lpts-1].c[0]/K0, KG[lpts-1].c[1]);
		fflush(outfile);

	}

	fclose(outfile);
}



// in-equilibrium derk
RETARRAY derk(RETARRAY (*func)(double l, double z[]), RETARRAY CalcKy, double Calcl, double Dl, unsigned n){
	int i,j;
	double Testdl, Testl[3];
	RETARRAY DKyDl[4], TestKy[3];

	DKyDl[0] = func(Calcl, CalcKy.c);
	for(i=0; i<3; i++){
		Testdl = Dl*((i+2)/2)/2;
		for(j=0; j<n; j++)
			TestKy[i].c[j] = CalcKy.c[j] + DKyDl[i].c[j]*Testdl;
		Testl[i] = Calcl + Dl*((i+1)/2)/2;
		DKyDl[i+1] = func(Testl[i], TestKy[i].c);
	}
	for(j=0; j<n; j++)
		CalcKy.c[j] += Dl*(DKyDl[0].c[j]+2.0*(DKyDl[1].c[j]+DKyDl[2].c[j])+DKyDl[3].c[j])/6.0;

	return CalcKy;
}


// out-of-equilibrium (Fokker Planck) derk
double derkK(double (*func)(double l, double z[]), RETARRAY CalcKG, double Calcl, double Dl){
	int i;
	double Testdl, Testl[3], DKDl[4];
	RETARRAY TestKG[3];

	DKDl[0] = func(Calcl, CalcKG.c);
	for(i=0; i<3; i++){
		Testdl = Dl*((i+2)/2)/2;
		TestKG[i].c[0] = CalcKG.c[0] + DKDl[i]*Testdl;
		TestKG[i].c[1] = CalcKG.c[1];
		Testl[i] = Calcl + Dl*((i+1)/2)/2;
		DKDl[i+1] = func(Testl[i], TestKG[i].c);
	}
	CalcKG.c[0] += Dl*(DKDl[0]+2*(DKDl[1]+DKDl[2])+DKDl[3])/6;

	return CalcKG.c[0];
}


// vortex-pair (KT) theory in-equilibrium recursion relations
RETARRAY EqRecRel(double l, double z[]){
	RETARRAY dzdl;

	dzdl.c[0] = -B*z[1]*z[1]*z[0]*z[0];
	dzdl.c[1] = (2-PI*z[0])*z[1];
	// dzdl.arr[2] = -2*PI*z[1]*exp(-2*l);

	return dzdl;
}


// vortex-pair (KT) theory out-of-equilibrium (Fokker Planck) recursion relation (for K)
double fpKRecRel(double l, double z[]){
	double dzdl;

	dzdl = -B*a04*z[0]*z[0]*exp(4*l)*z[1];

	return dzdl;
}


// vortex-pair (KT) theory out-of-equilibrium (Fokker Planck) recursion relation (for G)
double fpGRecRel(double l, double Dl, double zlo[], double zmd[], double zhi[]){
	double dzdt;

	//dzdt = exp(-2*l)*( PI*(zhi[0]*zhi[1]-zlo[0]*zlo[1])/Dl + (zhi[1]-2*zmd[1]+zlo[1])/(Dl*Dl) );
	dzdt = exp(-2*l)*( (zhi[0]*zhi[1]-zlo[0]*zlo[1])/2.0/Dl + (1/2.0/PI)*(zhi[1]-2*zmd[1]+zlo[1])/(Dl*Dl) );
	// Gstep[0] = exp(-2.0*x)*((cDiff/2.0/PI)*(Gin[j+1] - 2.0*Gin[j] + Gin[j-1])/(stepsize*stepsize) + (K[j+1]*G[j+1] - K[j-1]*G[j-1])/2.0/stepsize);

	return dzdt;
}



/*  Program Notes  */
/*=========================================================================*/
/*

See further documentation (UCLAresearchNotes.pdf, UCLAresearchPrograms.pdf, derk1.c) for elaboration.



*/