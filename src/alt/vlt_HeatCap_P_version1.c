/* filename: vlt_HeatCap_P_version1.c */
/* description: This C program calculates the molar specific heat capacity c_V of a system of liquid helium-4 at different temperatures, pressures, and length scales.
       The filename vlt_HeatCap.c refers to the "vortex loop theory" (vlt) of the superfluid phase transition and the molar specific "heat capacity" (HeatCap) of the system under examination.
 */
/* ext files: derk1.c */
/* author: Andrew Forrester <aforrester@ucla.edu> */

// NOTE: Be sure that N=3 in derk1.c when compiling.


/*  Function preparation  */
/*=========================================================================*/

#include <stdio.h>
#include <math.h>

#define IN
#define OUT
#define INOUT
#define PI     3.14159265358979323846
#define PISQ   9.86960440108935861883
#define A      4.0*PI*PI*PI/3.0
#define DK0    1e-07
#define P      0.71  // pressure in bars
#define N      3

// Parameters for Cc(P):
#define Cc0    1.10892106122653
#define Cc1   -0.0302090756207303
#define Cc2    0.000378628147396303
#define Cc3   -3.12142606957213e-06
#define Cc4    1.29173646324374e-08

// Parameters for K0c(P):
#define K0c0   0.295581776856052
#define K0c1   0.00411480553845432
#define K0c2   0.000209636865638176
#define K0c3  -6.16658606260075e-06
#define K0c4   9.20280978255458e-08

typedef struct {
	double l;
	double tempv;
	double thrmv[N];
} STATE;

typedef struct {
	double arr[N];
} RETARRAY;



/*  Function prototypes  */
/*=========================================================================*/
//   derk(func(),calculated,Dl,n)
//   ringrecrel(l,z)

extern STATE derk(RETARRAY (*func)(double l, double z[]), STATE calculated, double Dl, unsigned n);

RETARRAY ringrecrel(double l, double z[]);



/*  Function definitions  */
/*=========================================================================*/

main(){
	int i,j=0;
	double Cc,K0c,K0,Kr,cap,derivplus,derivminus;
	//double T,rho,rhos,x,Delok,a;
	STATE estimated, plus, minus;
	FILE *outfile;

	/* Set pressure dependence with the critical values of C and K0 */
	Cc  = 0.55; // Cc0  + Cc1*P  + Cc2*P*P  + Cc3*P*P*P  + Cc4*P*P*P*P;
	K0c = 4.646970466524e-01; // K0c0 + K0c1*P + K0c2*P*P + K0c3*P*P*P + K0c4*P*P*P*P;

	/* Print headings for data to screen and output file */
	outfile = fopen("vlt_HeatCap_P_24.80bar_DK0_1e-07_test.dat","w");  //  "vlt_HeatCap.out"
	fprintf(outfile,"%s\t%s\t%s\t%s\t%s\t\n", "tempv","cap","l","Kr/K0","F");
	printf(         "%s\t%s\t%s\t%s\t%s\t\n", "tempv","cap","l","Kr/K0","F");

	/* Loop through different temperatures (calculations for each temp are independent)? */
	for(i=0; i<35; i++){
		/* Decrement temperature (further below Tc) */
		estimated.tempv = (1.0e-8)*exp(0.4*i);

		/* Initialize lengthscale and energy values, zero in/decrements */
		estimated.l = 0.0;
		plus.l      = 0.0;
		minus.l     = 0.0;
		estimated.thrmv[2] = 0.0;
		plus.thrmv[2]      = 0.0;
		minus.thrmv[2]     = 0.0;

		/* Initialize temperature-dependent values at initial lengthscale */
		K0 = K0c/(1.0-estimated.tempv);
		estimated.thrmv[0] = K0;
		plus.thrmv[0]      = K0+DK0;
		minus.thrmv[0]     = K0-DK0;
		estimated.thrmv[1] = exp(-PISQ*K0*Cc);
		plus.thrmv[1]      = exp(-(K0+DK0)*PISQ*Cc);
		minus.thrmv[1]     = exp(-(K0-DK0)*PISQ*Cc);
		Kr = K0;

		while(1>0){
			estimated = derk(IN ringrecrel, INOUT estimated, IN 0.0001, IN N);
			plus      = derk(IN ringrecrel, INOUT plus,      IN 0.0001, IN N);
			minus     = derk(IN ringrecrel, INOUT minus,     IN 0.0001, IN N);

			Kr = estimated.thrmv[0]*exp(-estimated.l);
			j = j+1;
			derivplus = (plus.thrmv[2]-estimated.thrmv[2])/DK0;
			derivminus = (estimated.thrmv[2]-minus.thrmv[2])/DK0;
			cap = -8.31*K0*K0*(derivplus-derivminus)/DK0;
			/*if(j%2000==0)
				printf("%G %g %g %g \n%1.18e\n",t,Kr/K0,y[0],y[1],cap);
			if(j>10000)
				j=0;*/
			if(estimated.thrmv[0]>6)
				break;
		}
		cap = -8.31*K0*K0*(derivplus-derivminus)/DK0;

		fprintf(outfile,"%g\t%g\t%g\t%g\t%g\n",
			estimated.tempv, cap, estimated.l, Kr/K0, estimated.thrmv[2]);
		fflush(outfile);
		printf(         "%g\t%g\t%g\t%g\t%g\n",
			estimated.tempv, cap, estimated.l, Kr/K0, estimated.thrmv[2]);
	}

	fclose(outfile);
}



// vortex "ring" recursion relations
RETARRAY ringrecrel(double l, double z[]){
	RETARRAY dzdl;

	dzdl.arr[0] = z[0]-A*z[1]*z[0]*z[0];
	dzdl.arr[1] = z[1]*(6.0-PISQ*z[0]*(1.0-0.6*log(z[0])));
	dzdl.arr[2] = -PI*z[1]*exp(-3.0*l);

	return dzdl;
}


// deriv = de/dK0