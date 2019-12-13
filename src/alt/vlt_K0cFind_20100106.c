/*  File Comments  */
/*=========================================================================*/

/* filename: vlt_K0cFind.c */
/* description: This C program uses the vortex loop theory (vlt) of the superfluid phase transition to find the value of K0c, the critical value of K0, which is the bare dimensionless superfluidity ratio at the smallest length scale a0.
 */
/* ext files: derk1.c */
/* author: Andrew Forrester <aforrester@ucla.edu>  */

// NOTE: Be sure that N=2 in derk1.c when compiling.


/* RESULTS:
       Given  Cc = 1.03  and  K01 = 0.3,
       Previous calculated values have been K0c = 0.309146998265557
                                                  0.3091469984282580
                                                  0.3091469984216758
       The most updated calculation gets    K0c = 0.3091469984282580

Further Results:
Cc	K0c
1.20	0.2805443089652032
1.10	0.2964689389987730
1.06	0.3035407129346564
1.05	0.3053794236011743
1.04	0.3072479227697506
1.03	0.3091469984216724
1.02	0.3110774675059972
1.01	0.3130401773072950
1.00	0.3150360068926113
0.99	0.3170658686428607
0.98	0.3191307098748213
0.97	0.3212315145597378
0.90	0.3370353710663317
0.80	0.3636228655823348
0.70	0.3965669728967915
0.60	0.4387310072100324
0.55	0.4646970470439853
0.50	0.4951082899772941
0.40	0.5753434296585114
0.30	0.7010909420810506
0.20	0.9346549781449002
0.10	1.576187055119773
0.00	5.105045166053435

(See RhosAmps.ods for Cc(P), derived from Ahlers' 1973 data):
1.107411553625760	0.2952060040779711
1.060208920765020	0.3035026106143791
0.906690064087517	0.3354357910960570
0.730451367152917	0.3857194862980521
0.667516675997090	0.4090680770137908
0.588081513758998	0.4445731602870236
0.532586148563099	0.4747250509762295
*/



/*  Function Preparation  */
/*=========================================================================*/

#include <stdio.h>
#include <math.h>

#define IN
#define OUT
#define INOUT
#define PI        3.14159265358979323846
#define PISQ      9.86960440108935861883
#define A         4.0*PI*PI*PI/3.0
#define Cc        1.03
#define precision 1e-17
#define K01       0.3
// This value must be close to K0c for the program to find the actual K0c.
#define N         2

typedef struct {
	double l;
	double tempv;
	double thrmv[N];
} STATE;

typedef struct {
	double arr[N];
} RETARRAY;



/*  Function Prototypes  */
/*=========================================================================*/
//   derk(func(),calculated,Dl,n)
//   ringrecrel(l,z)

extern STATE derk(RETARRAY (*func)(double l, double z[]), STATE calculated, double Dl, unsigned n);

RETARRAY vltRecRel(double l, double z[]);



/*  Function Definitions  */
/*=========================================================================*/

main(){
	double K0,DK0,oldK,oldy;
	STATE estimated;

	/* Initial guess for K0 (should be close to K0c) and pre-initial adjustment of K0 */
	K0 = K01;
	DK0 = 0.1;  // It's a "pre-initial" adjustment because we change it before using it.

	/* Adjust K0 by finer and finer amounts until the adjustment DK0 reaches the desired precision */
	while(DK0>precision){
		/* Increase precision of adjustment of K0 */
		DK0 = DK0/2;

		/* Initialize lengthscale (and energy) values */
		estimated.l = 0.0;
		// estimated.thrmv[2] = 0.0;

		/* Initialize would-be temperature-dependent values at initial lengthscale */
		estimated.thrmv[0] = K0;
		estimated.thrmv[1] = exp(-PISQ*K0*Cc);
		oldK = estimated.thrmv[0];
		oldy = estimated.thrmv[1];

		/* Find whether K0 should be increased or decreased by DK0 to approach K0c */
		while(1==1){
			/* Calculate next estimated state, incrementing lengthscale */
			estimated = derk(IN vltRecRel, INOUT estimated, IN 0.0001, IN N);  // Dl=0.0001, n=N=2

			/* K is diverging, K0 too high */
			if(estimated.thrmv[0]>5.0) {
				//printf("K exploded\n");
				K0 = K0-DK0;
				break;
			}
			/* K is decreasing after l=6.0, K0 too low */
			if(estimated.l>6.0 && estimated.thrmv[0]<oldK) {
				//printf("K decreasing\n");
				K0 = K0+DK0;
				break;
			}

			/* If the above conditions have not yet been found, prepare for comparison with the next estimated state at a larger lengthscale */
			oldK = estimated.thrmv[0];
			oldy = estimated.thrmv[1];
		}
		printf("%1.15e %1.15e\n",K0,DK0);
	}
}



// vortex loop theory recursion relations
RETARRAY vltRecRel(double l, double z[]){
	RETARRAY dzdl;

	dzdl.arr[0] = z[0]-A*z[1]*z[0]*z[0];
	dzdl.arr[1] = z[1]*(6.0-PISQ*z[0]*(1.0-0.6*log(z[0])));
	// dzdl.arr[2] = -PI*z[1]*exp(-3.0*l);

	return dzdl;
}
