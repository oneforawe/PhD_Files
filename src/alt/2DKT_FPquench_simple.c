// 2DKT_FPquench_simple.c
#include <stdio.h>
#include <math.h>
#define PI       3.14159265358979323846
#define PISQ     9.86960440108935799230
#define PICU     31.00627668029981620634
#define B        4.0*PICU

/* Program parameters/inputs */
#define N        2
#define a0       1
#define a04      1  //  a0 to the fourth power
#define K0c      0.747853
#define dl       0.0001
#define lmax     1  //  10
// lpts = lmax/dl = 10/0.0001 = 100,000  (make sure this is divisible by ldatamax)
// the number of "l" points will really be lpts+1
#define ldatamax 100    //  max number of data points recorded per time step
// #define dt       1      //  time step (in units of "diffusion time"(?))
#define tmax     5      //  (later, 10000) max unitless time (see if(t==...) below)
#define iTfrac   0.95   //  initial temperature fraction Ti/Tkt
#define qTfrac   0.10   //  quench temperature fraction Tq/Tkt

typedef struct {
	double c[N];  //  c = "component/choose"
} RETARRAY;


/*  Function prototypes  */
/*=========================================================================*/
//  derk(func(),calcKy,calcl,Dl,n)
//  derkK(func(),calcKG,calcl,Dl,n)
//  derkG(func(),calcKG,calct,Dt,n);
//  EqRecRel(l,z)
//  FPRecRel(l,z)
RETARRAY derk(RETARRAY (*func)(double l, double z[]), RETARRAY calcKy, double calcl, double Dl, unsigned n);
double   derkK(double (*func)(double l, double z[]), RETARRAY calcKG, double calcl, double Dl);
//double*  derkG(double (*func)(double l, double Dl, double zlo[], double zmd[], double zhi[]), RETARRAY calcKG[], double calcl[], double Dl, double Dt);
RETARRAY EqRecRel(double l, double z[]);
double   fpKRecRel(double l, double z[]);
double   fpGRecRel(double l, double Dl, double zlo[], double zmd[], double zhi[]);
//double   derkG(double (*func)(double l, double Dl, double zlo[], double zmd[], double zhi[]), RETARRAY calcKGlo, RETARRAY calcKGmd, RETARRAY calcKGhi, double calcl, double Dl, double Dt);


/*  Function definitions  */
/*=========================================================================*/

main(){
	/* Main function definitions */
	int t,i,lsteps,lpts;
	double tau,K0,y0;
	lsteps = lmax/dl;
	lpts = lsteps+1;  // from l=0 to l=lmax, inclusive
	double l[lpts], dGdt[lpts];
	RETARRAY Ky[lpts], KG[lpts+1];
	KG[lpts].c[0] = KG[lpts].c[1] = 0;  // at l=lmax+dl, enforcing boundary condition (correct way of saying it?)
	FILE *outfile;
	char filename[100];

	int j;
	double testDGDt[lpts+1][4], testdt;
	double k1[lpts+1], k2[lpts+1], k3[lpts+1], k4[lpts+1];
	RETARRAY testKG[lpts+1][3];

	double dt=1,dtLimit;
//	double finalDGDt[lpts];

	/* Print headings for data to screen and output file */
	sprintf(filename, "2DKT_FPquench_simple_lmax%i_K0c%7.6f.dat", lmax, K0c);
	outfile = fopen(filename,"w");  //  E.g., "2DKT_FPquench.out" or "2DKT_FPquench_lmax100_K0c0.747853.dat"
	fprintf(outfile,"%s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","G","K/K0");
	printf(         "%s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","G","K/K0");

    /*** Initial (t=0), Equilibrium Calculations ***/

	/* Initialize */
	tau = 1.0-iTfrac;
	l[0] = 0.0;
	K0 = K0c/(1.0-tau);
	y0 = exp(-PISQ*K0/2);
	Ky[0].c[0] = K0;  //  K
	KG[0].c[0] = K0;  //  K
	Ky[0].c[1] = y0;  //  y
	KG[0].c[1] = exp(-4*l[0])*y0*y0/a04;  // G (not = 0.0)

	/* Print initial data to screen and output file */
	fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\n", 1-tau, tau, l[0], KG[0].c[0], KG[0].c[1], KG[0].c[0]/K0);
	printf(         "%g\t%g\t%g\t%g\t%g\t%g\n", 1-tau, tau, l[0], KG[0].c[0], KG[0].c[1], KG[0].c[0]/K0);

	/* Calculate state-variables at progressively larger length scales, up to coherence length, via Runge-Kutta method (derk) */
	for(i=1; i<lpts; i++){
		/* Calculate next Ky and KG state, incrementing length scale */
		Ky[i] = derk(EqRecRel, Ky[i-1], l[i-1], dl, N);
		l[i] = l[i-1] + dl;
		KG[i].c[0] = Ky[i].c[0];
		KG[i].c[1] = exp(-4*l[i])*Ky[i].c[1]*Ky[i].c[1]/a04;

		/* Print at most [ldatamax] (e.g., 1000) data points to screen and output file */
		if(i%(lsteps/ldatamax)==0||i<10){
			fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\n", 1-tau, tau, l[i], KG[i].c[0], KG[i].c[1], KG[i].c[0]/K0);
			printf(         "%g\t%g\t%g\t%g\t%g\t%g\n", 1-tau, tau, l[i], KG[i].c[0], KG[i].c[1], KG[i].c[0]/K0);
		}
	}

    /*** After Temperature Drop/Quench ***/

	/* Define the timestep size */
	dtLimit = 2.5*(dl*dl)/sqrt((4.0/PI/PI) + K0*dl*K0*dl);
	dt = 0.9*dtLimit;

	/* Advance in time */
	for(t=1; t<=tmax; t++){

		/* Initialize */
		tau = 1.0-qTfrac;
		// l[0] = 0.0;  //  (already done)
		K0 = K0c/(1.0-tau);
		y0 = exp(-PISQ*K0/2);
		KG[0].c[0] = K0;  //  K
		KG[0].c[1] = exp(-4*l[0])*y0*y0/a04;  // G (not = 0.0)  //  (correct location?)

		/* Calculate K at progressively larger length scales, up to coherence length, via Runge-Kutta method (derk1) */
		for(i=1; i<lpts; i++)
			KG[i].c[0] = derkK(fpKRecRel, KG[i-1], l[i-1], dl);

		/* Calculate dGdt, for all length scales */
/*		dGdt = derkG(fpGRecRel, KG, l, dl, dt);       */
/*  */		for(j=0; j<4; j++){
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
			dGdt[i] = (testDGDt[i][0]+2*(testDGDt[i][1]+testDGDt[i][2])+testDGDt[i][3])/6;

		/* Step G in time, for all length scales */
		for(i=1; i<lpts; i++)
			KG[i].c[1] += dGdt[i]*dt;


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
		//testKG[0][0].c[1] = testKG[1][0].c[1];
		//testKG[lpts-1][0].c[1] = testKG[lpts-2][0].c[1];

		for(i=1; i<lpts-1; i++)
			testDGDt[i][1]    = fpGRecRel(l[i], dl, testKG[i-1][0].c, testKG[i][0].c, testKG[i+1][0].c);
		for(i=0; i<lpts; i++){
			testKG[i][1].c[0] = KG[i].c[0];
			k2[i]             = testDGDt[i][1]*dt;
			testKG[i][1].c[1] = KG[i].c[1] - k1[i] + k2[i];
		}
		//testKG[0][1].c[1] = testKG[1][1].c[1];
		//testKG[lpts-1][1].c[1] = testKG[lpts-2][1].c[1];

		for(i=1; i<lpts-1; i++)
			testDGDt[i][2]    = fpGRecRel(l[i], dl, testKG[i-1][1].c, testKG[i][1].c, testKG[i+1][1].c);
		for(i=0; i<lpts; i++){
			testKG[i][2].c[0] = KG[i].c[0];
			k3[i]             = testDGDt[i][2]*dt;
			testKG[i][2].c[1] = KG[i].c[1] + k1[i]/6.0 + 2.0*k2[i]/3.0 + k3[i]/6.0;
		}
		//testKG[0][2].c[1] = testKG[1][2].c[1];
		//testKG[lpts-1][2].c[1] = testKG[lpts-2][2].c[1];

		for(i=1; i<lpts-1; i++)
			testDGDt[i][3]    = fpGRecRel(l[i], dl, testKG[i-1][2].c, testKG[i][2].c, testKG[i+1][2].c);
		for(i=0; i<lpts; i++){
			k4[i]             = testDGDt[i][3]*dt;
			KG[i].c[1]        = KG[i].c[1] + k1[i]/6.0 + k2[i]/3.0 + k3[i]/3.0 + k4[i]/6.0;
		}  */


		/* At certain times, print out the results */
		if(t<5){   //  t==1||t==10||t==100||t==1000||t==10000){  //  (log10(t)%1==0){
			fprintf(outfile,"\n%s%i\n\n", "t = ",t);
			fprintf(outfile,"%s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","G","K/K0");
			for(i=0; i<lpts; i++){
				/* Print at most [ldatamax] (e.g., 1000) data points to screen and output file */
				if(i%(lsteps/ldatamax)==0||i<10){
					fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\n", 1-tau, tau, l[i], KG[i].c[0], KG[i].c[1], KG[i].c[0]/K0);
					printf(         "%g\t%g\t%g\t%g\t%g\t%g\n", 1-tau, tau, l[i], KG[i].c[0], KG[i].c[1], KG[i].c[0]/K0);
					fflush(outfile);
				}
			}
		}
	}

	fclose(outfile);
}



// in-equilibrium derk
RETARRAY derk(RETARRAY (*func)(double l, double z[]), RETARRAY calcKy, double calcl, double Dl, unsigned n){
	int i,j;
	double testdl, testl[3];
	RETARRAY DKyDl[4], testKy[3];

	DKyDl[0] = func(calcl, calcKy.c);
	for(i=0; i<3; i++){
		testdl = Dl*((i+2)/2)/2;
		for(j=0; j<n; j++)
			testKy[i].c[j] = calcKy.c[j] + DKyDl[i].c[j]*testdl;
		testl[i] = calcl + Dl*((i+1)/2)/2;
		DKyDl[i+1] = func(testl[i], testKy[i].c);
	}
	for(j=0; j<n; j++)
		calcKy.c[j] += Dl*(DKyDl[0].c[j]+2*(DKyDl[1].c[j]+DKyDl[2].c[j])+DKyDl[3].c[j])/6;

	return calcKy;
}


// out-of-equilibrium (Fokker Planck) derk
double derkK(double (*func)(double l, double z[]), RETARRAY calcKG, double calcl, double Dl){
	int i;
	double testdl, testl[3], DKDl[4];
	RETARRAY testKG[3];

	DKDl[0] = func(calcl, calcKG.c);
	for(i=0; i<3; i++){
		testdl = Dl*((i+2)/2)/2;
		testKG[i].c[0] = calcKG.c[0] + DKDl[i]*testdl;
		testKG[i].c[1] = calcKG.c[1];
		testl[i] = calcl + Dl*((i+1)/2)/2;
		DKDl[i+1] = func(testl[i], testKG[i].c);
	}
	calcKG.c[0] += Dl*(DKDl[0]+2*(DKDl[1]+DKDl[2])+DKDl[3])/6;

	return calcKG.c[0];
}


// out-of-equilibrium (Fokker Planck) derk
/* double (*derkG)(double (*func)(double l, double Dl, double zlo[], double zmd[], double zhi[]), RETARRAY calcKG[], double calcl[], double Dl, double Dt){
	int j,i;
	double testDGDt[lpts][4], finalDGDt[lpts], testdt;
	RETARRAY testKG[lpts][3];


	for(j=1; j<lpts-1; j++)
		testDGDt[j][0] = fpGRecRel(calcl[j], Dl, calcKG[j-1].c, calcKG[j].c, calcKG[j+1].c);
	for(i=0; i<3; i++){
		testdt = Dt*((i+2)/2)/2;
		for(j=1; j<lpts-1; j++){
			testKG[j][i].c[0] = calcKG[j].c[0];
			testKG[j][i].c[1] = calcKG[j].c[1] + testDGDt[j][i]*testdt;
		}
		for(j=1; j<lpts-1; j++)
			testDGDt[j][i+1] = fpGRecRel(calcl[j], Dl, testKG[j-1][i].c, testKG[j][i].c, testKG[j+1][i].c);
	}
	for(j=1; j<lpts-1; j++)
		finalDGDt[j] = (testDGDt[j][0]+2*(testDGDt[j][1]+testDGDt[j][2])+testDGDt[j][3])/6;

	return finalDGDt;
} */


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

	dzdt = exp(-2*l)*( PI*(zhi[0]*zhi[1]-zlo[0]*zlo[1])/Dl + (zhi[1]-2*zmd[1]+zlo[1])/(Dl*Dl) );
	// dzdt = exp(-2*l)*( (zhi[0]*zhi[1]-zlo[0]*zlo[1])/2.0/Dl + (1/2.0/PI)*(zhi[1]-2*zmd[1]+zlo[1])/(Dl*Dl) );
	// Gstep[0] = exp(-2.0*x)*((cDiff/2.0/PI)*(Gin[j+1] - 2.0*Gin[j] + Gin[j-1])/(stepsize*stepsize) + (K[j+1]*G[j+1] - K[j-1]*G[j-1])/2.0/stepsize);

	return dzdt;
}












		/* Calculate dGdt, for all length scales */  // Simple derivative time-stepping doesn't work.  Need to use derk method.
		//for(i=1; i<lpts-1; i++)
		//	dGdt[i] = exp(-2*l[i])*( PI*(KG[i+1].c[0]*KG[i+1].c[1]-KG[i-1].c[0]*KG[i-1].c[1])/dl + (KG[i+1].c[1]-2*KG[i].c[1]+KG[i-1].c[1])/(dl*dl) );
		//dGdt[lpts-1] = exp(-2*l[lpts-1])*( PI*(-KG[lpts-2].c[0]*KG[lpts-2].c[1])/dl + (-2*KG[lpts-1].c[1]+KG[lpts-2].c[1])/(dl*dl) );

		/* For bug detection, to check that calculations make sense...
		fprintf(outfile,"\n\n%s\t%s\t%s\n", "l[]", "K[]", "G[]");
		for(i=0; i<10; i++){
			fprintf(outfile,"%g\t%g\t%g\n", l[i], KG[i].c[0], KG[i].c[1]);
			fflush(outfile);
		}
		for(i=1; i<10; i++){
			fprintf(outfile,"\n%s\n", "dGdt[i] = exp(-2*l[i])*( PI*( K[i+1]*G[i+1] - K[i-1]*G[i-1] )/dl                           + ( G[i+1] - 2*G[i] + G[i-1] )/(dl*dl) )");
			fprintf(outfile,"%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s%g%s\n", "dGdt[i] = exp(-2*{",l[i],"})*( PI*({",KG[i+1].c[0],"}*{",KG[i+1].c[1],"}-{",KG[i-1].c[0],"}*{",KG[i-1].c[1],"})/{",dl,"} + ({",KG[i+1].c[1],"}-2*{",KG[i].c[1],"}+{",KG[i-1].c[1],"})/({",dl,"}*{",dl,"}) )");
			fprintf(outfile,"%s%g%s%g%s%g%s%g%s%g%s\n", "dGdt[i] = {",exp(-2*l[i]),"*( {", PI*(KG[i+1].c[0]*KG[i+1].c[1]-KG[i-1].c[0]*KG[i-1].c[1]),"}/{",dl,"} + {",(KG[i+1].c[1]-2*KG[i].c[1]+KG[i-1].c[1]),"}/{",dl*dl,"} )");
			fprintf(outfile,"%s%g%s%g%s%g%s\n", "dGdt[i] = {",exp(-2*l[i]),"*( {", PI*(KG[i+1].c[0]*KG[i+1].c[1]-KG[i-1].c[0]*KG[i-1].c[1])/dl,"} + {",(KG[i+1].c[1]-2*KG[i].c[1]+KG[i-1].c[1])/(dl*dl),"} )");
			fprintf(outfile,"%s%g\n", "dGdt[i] = {",exp(-2*l[i])*( PI*(KG[i+1].c[0]*KG[i+1].c[1]-KG[i-1].c[0]*KG[i-1].c[1])/dl + (KG[i+1].c[1]-2*KG[i].c[1]+KG[i-1].c[1])/(dl*dl) ) );
			fprintf(outfile,"%s%g\n", "dGdt[i] = ", dGdt[i]);
			fflush(outfile);
		} */



// out-of-equilibrium (Fokker Planck) derk
/* Wrong order for derk procedure
double derkG(double (*func)(double l, double Dl, double zlo[], double zmd[], double zhi[]), RETARRAY calcKGlo, RETARRAY calcKGmd, RETARRAY calcKGhi, double calcl, double Dl, double Dt){
	int i;
	double testdt, testt[3], DGDt[4], DGDtFinal;
	RETARRAY testKGlo[3], testKGmd[3], testKGhi[3];

	DGDt[0] = func(calcl, Dl, calcKGlo.c, calcKGmd.c, calcKGhi.c);
	for(i=0; i<3; i++){
		testKGlo[i].c[0] = calcKGlo.c[0];
		testKGmd[i].c[0] = calcKGmd.c[0];
		testKGhi[i].c[0] = calcKGhi.c[0];
		testdt = Dt*((i+2)/2)/2;
		testKGlo[i].c[1] = calcKGlo.c[1] + DGDt[i]*testdt;
		testKGmd[i].c[1] = calcKGmd.c[1] + DGDt[i]*testdt;
		testKGhi[i].c[1] = calcKGhi.c[1] + DGDt[i]*testdt;
		DGDt[i+1] = func(calcl, Dl, testKGlo[i].c, testKGmd[i].c, testKGhi[i].c);
	}
	DGDtFinal = (DGDt[0]+2*(DGDt[1]+DGDt[2])+DGDt[3])/6;

	return DGDtFinal;
} */




/*

	for(j=1; j<lpts-1; j++)
		testDGDt[j][0] = fpGRecRel(l[j], dl, KG[j-1].c, KG[j].c, KG[j+1].c);
	for(i=0; i<3; i++){
		testdt = dt*((i+2)/2)/2;
		for(j=1; j<lpts-1; j++){
			testKG[j][i].c[0] = KG[j].c[0];
			testKG[j][i].c[1] = KG[j].c[1] + testDGDt[j][i]*testdt;
		}
		for(j=1; j<lpts-1; j++)
			testDGDt[j][i+1] = fpGRecRel(l[j], dl, testKG[j-1][i].c, testKG[j][i].c, testKG[j+1][i].c);
	}
	for(j=1; j<lpts-1; j++)
		finalDGDt[j] = (testDGDt[j][0]+2*(testDGDt[j][1]+testDGDt[j][2])+testDGDt[j][3])/6;


*/