/*  File Comments  *//*=========================================================================*//* FILENAME: vlt_HeatCap_P.c *//* AUTHOR: Andrew Forrester <aforrester@ucla.edu> *//* DESCRIPTION:   * This C program calculates the (length-scale-dependent) molar specific heat capacity at constant pressure c_P (cp) (and other thermodynamic quantities) of a system of liquid helium-4 (4He) at various temperatures and length scales, given a pressure (in the parameters/inputs below).   * The system is a spherical region (of variable diameter greater than D) of liquid 4He at some constant pressure, near or below the critical transition temperature Tc (which depends on pressure) separating normal fluid and superfluid (aka the lambda temperature).   * The program uses the vortex-loop theory of the superfluid phase transition (i.e., its differential recursion relations) of Gary A. Williams and Subodh R. Shenoy and a 4th-order Runge-Kutta integration method to calculate (length-scale-dependent) thermodynamic quantities (especially cp) of the system at a set of temperatures {Tn}, at various length scales (no greater than a length scale corresponding to the the limit K<Kthresh), and at some pressure (which also determines the choice of Cc, the critical dimensionless vortex-loop-core energy, and K0c, the critical "bare" dimensionless superfluid ratio).   * In the basename "vlt_HeatCap_P", "vlt" refers to Vortex Loop Theory (of the superfluid phase transition), "HeatCap" refers to and the molar specific HEAT CAPacity of the system under examination, and "P" refers to the pressure, which must be specified in units of bars in the parameters/inputs below.   * This program could be subsumed in vlt_ThermStates.c...   Inputs:   * Kthresh -   * Cc      -   * K0c     -   * P       -   * dl      -   * DK0     - This allows you to take derivatives with respect to K0 (corresponds to derivatives wrt temperature).   * digits  - This allows you to decide when "eimprecise" is subtracted from "e", leaving "de" for the remaining calculations.  (Since a double float numbers only have about 16 digits, we have to use a trick to calculate more digits in e and make the derivative of e with respect to K0 more precise.  In the "First pass", we calculate eimprecise, which is e to about 15 digits.  In the "Second pass", as the calculated value of e stabilizes, eimprecise is subtraced off of e, leaving what we call "de".  Since only differences in e matter, we continue the calculations with de.  The value of de will inevitably be on the order of 10^-15 in the end, but digits influences when the de values start to be used and thus has some influence on the precision of the calculations.)   * Dexp    - Why not?   * You could add options for selecting different sets of temperatures.  Right now, you can just "Decrement temperature (further below Tc)".   Output:   * This program returns text file of data in a two dimensional array.*//* EXT FILES: vlt_derk.c *//* COMPILE NOTES:   * To run, type "gcc -lm vlt_HeatCap_P.c vlt_derk.c" without the quotes.   * NOTE: Be sure that N=3 in vlt_derk.c when compiling.   * NOTE: Update P, Cc, K0c, and the output filename before running.*//*  Function Preparation  *//*=========================================================================*//* Standard routine header files */#include <stdio.h>#include <math.h>/* Constants and labels */#define IN      // input label#define OUT     // output label#define INOUT   // input/output label#define PI      3.14159265358979323846#define PISQ    9.86960440108935861883#define A       4.0*PI*PI*PI/3.0#define R       8.31447215#define N       3/* Program parameters/inputs */#define Kthresh 5.0#define P       0.050  // pressure in bars#define Cc      1.105023144261770#define K0c     0.295611433115442#define dl      0.0001#define DK0     1e-08#define digits  8//#define Dexp/*(See RhosAmps.ods for Cc(P), derived from Ahlers' 1973 data):P (bar)	Cc(P)	K0c(Cc) (from K0cFind.c)0.050	1.105023144261770	0.2956114331154421.646	1.058195551375160	0.3038703306608067.328	0.905924430200288	0.33561783603736415.031	0.730725772346878	0.38562526579595118.180	0.667963030000677	0.40888923739440222.533	0.588565898089815	0.44433174497922325.868	0.532988906023542	0.474486704164683#define P       0.050  // pressure in bars#define Cc      1.105023144261770#define K0c     0.295611433115442#define P       1.646  // pressure in bars#define Cc      1.058195551375160#define K0c     0.303870330660806#define P       7.328  // pressure in bars#define Cc      0.905924430200288#define K0c     0.335617836037364#define P       15.031  // pressure in bars#define Cc      0.730725772346878#define K0c     0.385625265795951#define P       18.180  // pressure in bars#define Cc      0.667963030000677#define K0c     0.408889237394402#define P       22.533  // pressure in bars#define Cc      0.588565898089815#define K0c     0.444331744979223#define P       25.868  // pressure in bars#define Cc      0.532988906023542#define K0c     0.474486704164683*//* // Cc(P) fit parameters:#define Cc00    1.10652102266763#define Cc01   -0.0299766068265842#define Cc02    0.000380954449384652#define Cc03   -3.60465470172711e-06#define Cc04   -1.13928690475838e-08#define Cc05    5.27191245296858e-09#define Cc06   -4.02970077711438e-10#define Cc07    1.8784597083754e-11#define Cc08   -5.69181439968371e-13#define Cc09    1.12069994170299e-14#define Cc10   -1.38157107961357e-16#define Cc11    9.67488197927817e-19#define Cc12   -2.93300475897348e-21// (See RhosAmps.ods)// K0c(P) fit parameters:// USE vlt_K0cFind.c instead!#define K0c00   0.295357018625725#define K0c01   0.00508602480223517#define K0c02   5.14719279164173e-05#define K0c03   4.7417125185559e-07#define K0c04   2.15616844605396e-08#define K0c05  -1.76018997011296e-09#define K0c06   1.31251949946304e-10#define K0c07  -6.00036335403098e-12#define K0c08   1.82312804922308e-13#define K0c09  -3.62344372128926e-15#define K0c10   4.56051595849655e-17#define K0c11  -3.30108155723952e-19#define K0c12   1.05873049007624e-21// (See RhosAmps.ods)*/// Parameters for Tc(P):#define Tc0    2.173494256#define Tc1   -0.009824996#define Tc2   -0.000118194#define Tc3   -0.000000437#define Tc4    0.000000007// Parameters for rho(T,P):#define rho00    145.145109496329000#define rho10   -0.097653969305915#define rho20    0.334163407001684#define rho30   -0.446930785976304#define rho40    0.181879478545246#define rho01    1.744776044955830#define rho11   -0.091953899317905#define rho21    0.179844560873926#define rho31   -0.133606331352667#define rho41    0.041022551424992#define rho02   -0.049165537969017#define rho12    0.007106988980704#define rho22   -0.008230542254959#define rho32    0.000609542602247#define rho42    0.001149167753923#define rho03    0.001341503764375#define rho13   -0.000362007479156#define rho23    0.000358809384119#define rho33    0.000064818395436#define rho43   -0.000104112551303#define rho04   -0.000016990729415#define rho14    0.000005538203683#define rho24   -0.000003157734111#define rho34   -0.000004999673069#define rho44    0.000003413312235// Parameters for alpha(T,P):#define alpha00   -0.00132545289717059#define alpha10   0.0179528212646871#define alpha20   -0.077814417132819#define alpha30   0.148812649901035#define alpha40   -0.135348183536229#define alpha50   0.0575865394848149#define alpha60   -0.00942356361818271#define alpha01   0.00169244447357293#define alpha11   -0.021665108348567#define alpha21   0.0875997904161722#define alpha31   -0.155075446681196#define alpha41   0.133241381828243#define alpha51   -0.0547135426838175#define alpha61   0.00857815766443886#define alpha02   -0.000402016588457985#define alpha12   0.00500646576912193#define alpha22   -0.0196997735925578#define alpha32   0.0343311248462808#define alpha42   -0.0295791060596021#define alpha52   0.0123277022988963#define alpha62   -0.00198704719059122#define alpha03   3.50822131497139e-05#define alpha13   -0.000417434943625329#define alpha23   0.00157100359153513#define alpha33   -0.0026926366274236#define alpha43   0.00238134017220718#define alpha53   -0.00104902823471905#define alpha63   0.000181757286633977#define alpha04   -1.51074156982117e-06#define alpha14   1.69620331490882e-05#define alpha24   -6.10050371264495e-05#define alpha34   0.000107369015843627#define alpha44   -0.000106541517748163#define alpha54   5.44104698460238e-05#define alpha64   -1.08691954908159e-05#define alpha05   3.11237240810798e-08#define alpha15   -3.20093978259846e-07#define alpha25   1.08841581501832e-06#define alpha35   -2.12545135194623e-06#define alpha45   2.63498247456092e-06#define alpha55   -1.62910668019171e-06#define alpha65   3.71525294733195e-07#define alpha06   -2.37222132465804e-10#define alpha16   2.08583118201288e-09#define alpha26   -6.5282210389408e-09#define alpha36   1.6933030216303e-08#define alpha46   -2.89786898905826e-08#define alpha56   2.12469804642866e-08#define alpha66   -5.29362586788296e-09/* Data type definitions */typedef struct {	double l;	double tempv;	double thrmv[N];} STATE;typedef struct {	double arr[N];} RETARRAY;/*  Function Prototypes  *//*=========================================================================*///   derk(func(),calculated,Dl,n)//   vltRecRel(l,z)extern STATE derk(RETARRAY (*func)(double l, double z[]), STATE calculated, double Dl, unsigned n);RETARRAY vltRecRel(double l, double z[]);/*  Function Definitions  *//*=========================================================================*/main(){	int i,j=0;	double K0,Kr,Tc,T,rho,DrhoDT,D2rhoDT2,a0,Da0DT,D2a0DT2;	double DeDK0plus,DeDK0minus,DeDK0,D2eDK02,alpha,DalphaDT,cp;	double eimprecise;	//double Cc,K0c,rhos,x,Delok,a;	STATE estimated, plus, minus;	FILE *outfile;	/* Set pressure dependence with the critical values of C and K0 */	// Cc  = Cc00 + Cc01*P + Cc02*P*P + Cc03*P*P*P + Cc04*P*P*P*P + Cc05*P*P*P*P*P + Cc06*P*P*P*P*P*P + Cc07*P*P*P*P*P*P*P + Cc08*P*P*P*P*P*P*P*P + Cc09*P*P*P*P*P*P*P*P*P + Cc10*P*P*P*P*P*P*P*P*P*P + Cc11*P*P*P*P*P*P*P*P*P*P*P + Cc12*P*P*P*P*P*P*P*P*P*P*P*P;	// K0c = K0c00 + K0c01*P + K0c02*P*P + K0c03*P*P*P + K0c04*P*P*P*P + K0c05*P*P*P*P*P + K0c06*P*P*P*P*P*P + K0c07*P*P*P*P*P*P*P + K0c08*P*P*P*P*P*P*P*P + K0c09*P*P*P*P*P*P*P*P*P + K0c10*P*P*P*P*P*P*P*P*P*P + K0c11*P*P*P*P*P*P*P*P*P*P*P + K0c12*P*P*P*P*P*P*P*P*P*P*P*P;	/* Open/name output file, and print headings for data to screen and output file */	outfile = fopen("vlt_HeatCap_P_00.050_DK0_5e-08.dat","w");  //  E.g., "vlt_HeatCap_P.out" or "vlt_HeatCap_P_00.050_DK0_1e-07.dat"	fprintf(outfile,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",		"T/Tc","1-T/Tc","l","K","y","e","Kr","Kr/K0","cp");	printf(         "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",		"T/Tc","1-T/Tc","l","K","y","e","Kr","Kr/K0","cp");	/* Loop through different temperatures (calculations for each temp are independent) */	for(i=0; i<35; i++){		/* Decrement temperature (further below Tc) */		estimated.tempv = (1.0e-8)*exp(0.4*i);		/* Set temperature- and pressure-dependent values */		Tc = Tc0 + Tc1*P + Tc2*P*P + Tc3*P*P*P + Tc4*P*P*P*P;		T = Tc*(1-estimated.tempv);		rho = rho00         + rho10*T         + rho20*T*T         + rho30*T*T*T         + rho40*T*T*T*T 		    + rho01*P       + rho11*T*P       + rho21*T*T*P       + rho31*T*T*T*P       + rho41*T*T*T*T*P 		    + rho02*P*P     + rho12*T*P*P     + rho22*T*T*P*P     + rho32*T*T*T*P*P     + rho42*T*T*T*T*P*P 		    + rho03*P*P*P   + rho13*T*P*P*P   + rho23*T*T*P*P*P   + rho33*T*T*T*P*P*P   + rho43*T*T*T*T*P*P*P 		    + rho04*P*P*P*P + rho14*T*P*P*P*P + rho24*T*T*P*P*P*P + rho34*T*T*T*P*P*P*P + rho44*T*T*T*T*P*P*P*P;		DrhoDT = rho10         + 2*rho20*T         + 3*rho30*T*T         + 4*rho40*T*T*T 		       + rho11*P       + 2*rho21*T*P       + 3*rho31*T*T*P       + 4*rho41*T*T*T*P 		       + rho12*P*P     + 2*rho22*T*P*P     + 3*rho32*T*T*P*P     + 4*rho42*T*T*T*P*P 		       + rho13*P*P*P   + 2*rho23*T*P*P*P   + 3*rho33*T*T*P*P*P   + 4*rho43*T*T*T*P*P*P 		       + rho14*P*P*P*P + 2*rho24*T*P*P*P*P + 3*rho34*T*T*P*P*P*P + 4*rho44*T*T*T*P*P*P*P;		D2rhoDT2 = 2*rho20         + 6*rho30*T         + 12*rho40*T*T 		         + 2*rho21*P       + 6*rho31*T*P       + 12*rho41*T*T*P 		         + 2*rho22*P*P     + 6*rho32*T*P*P     + 12*rho42*T*T*P*P 		         + 2*rho23*P*P*P   + 6*rho33*T*P*P*P   + 12*rho43*T*T*P*P*P 		         + 2*rho24*P*P*P*P + 6*rho34*T*P*P*P*P + 12*rho44*T*T*P*P*P*P;		a0 = 549*Tc*K0c/rho;  // in angstroms		// a0(T,P) = ( m^2 kB Tc(P) K0c(P) ) / (hbar^2 rho(T,P) )		// where m = 6.65e-27 is the mass of 4He atom in kg, kB is the Boltzmann constant, and hbar is the Dirac constant		Da0DT = -549*Tc*K0c*DrhoDT/(rho*rho);		D2a0DT2 = 2*549*Tc*K0c*DrhoDT*DrhoDT/(rho*rho*rho) - 549*Tc*K0c*D2rhoDT2/(rho*rho);		alpha = alpha00             + alpha10*T             + alpha20*T*T             + alpha30*T*T*T             + alpha40*T*T*T*T             + alpha50*T*T*T*T*T             + alpha60*T*T*T*T*T*T 		      + alpha01*P           + alpha11*T*P           + alpha21*T*T*P           + alpha31*T*T*T*P           + alpha41*T*T*T*T*P           + alpha51*T*T*T*T*T*P           + alpha61*T*T*T*T*T*T*P 		      + alpha02*P*P         + alpha12*T*P*P         + alpha22*T*T*P*P         + alpha32*T*T*T*P*P         + alpha42*T*T*T*T*P*P         + alpha52*T*T*T*T*T*P*P         + alpha62*T*T*T*T*T*T*P*P 		      + alpha03*P*P*P       + alpha13*T*P*P*P       + alpha23*T*T*P*P*P       + alpha33*T*T*T*P*P*P       + alpha43*T*T*T*T*P*P*P       + alpha53*T*T*T*T*T*P*P*P       + alpha63*T*T*T*T*T*T*P*P*P 		      + alpha04*P*P*P*P     + alpha14*T*P*P*P*P     + alpha24*T*T*P*P*P*P     + alpha34*T*T*T*P*P*P*P     + alpha44*T*T*T*T*P*P*P*P     + alpha54*T*T*T*T*T*P*P*P*P     + alpha64*T*T*T*T*T*T*P*P*P*P 		      + alpha05*P*P*P*P*P   + alpha15*T*P*P*P*P*P   + alpha25*T*T*P*P*P*P*P   + alpha35*T*T*T*P*P*P*P*P   + alpha45*T*T*T*T*P*P*P*P*P   + alpha55*T*T*T*T*T*P*P*P*P*P   + alpha65*T*T*T*T*T*T*P*P*P*P*P 		      + alpha06*P*P*P*P*P*P + alpha16*T*P*P*P*P*P*P + alpha26*T*T*P*P*P*P*P*P + alpha36*T*T*T*P*P*P*P*P*P + alpha46*T*T*T*T*P*P*P*P*P*P + alpha56*T*T*T*T*T*P*P*P*P*P*P + alpha66*T*T*T*T*T*T*P*P*P*P*P*P;		DalphaDT = alpha10             + 2*alpha20*T             + 3*alpha30*T*T             + 4*alpha40*T*T*T             + 5*alpha50*T*T*T*T             + 6*alpha60*T*T*T*T*T 		         + alpha11*P           + 2*alpha21*T*P           + 3*alpha31*T*T*P           + 4*alpha41*T*T*T*P           + 5*alpha51*T*T*T*T*P           + 6*alpha61*T*T*T*T*T*P 		         + alpha12*P*P         + 2*alpha22*T*P*P         + 3*alpha32*T*T*P*P         + 4*alpha42*T*T*T*P*P         + 5*alpha52*T*T*T*T*P*P         + 6*alpha62*T*T*T*T*T*P*P 		         + alpha13*P*P*P       + 2*alpha23*T*P*P*P       + 3*alpha33*T*T*P*P*P       + 4*alpha43*T*T*T*P*P*P       + 5*alpha53*T*T*T*T*P*P*P       + 6*alpha63*T*T*T*T*T*P*P*P 		         + alpha14*P*P*P*P     + 2*alpha24*T*P*P*P*P     + 3*alpha34*T*T*P*P*P*P     + 4*alpha44*T*T*T*P*P*P*P     + 5*alpha54*T*T*T*T*P*P*P*P     + 6*alpha64*T*T*T*T*T*P*P*P*P 		         + alpha15*P*P*P*P*P   + 2*alpha25*T*P*P*P*P*P   + 3*alpha35*T*T*P*P*P*P*P   + 4*alpha45*T*T*T*P*P*P*P*P   + 5*alpha55*T*T*T*T*P*P*P*P*P   + 6*alpha65*T*T*T*T*T*P*P*P*P*P 		         + alpha16*P*P*P*P*P*P + 2*alpha26*T*P*P*P*P*P*P + 3*alpha36*T*T*P*P*P*P*P*P + 4*alpha46*T*T*T*P*P*P*P*P*P + 5*alpha56*T*T*T*T*P*P*P*P*P*P + 6*alpha66*T*T*T*T*T*P*P*P*P*P*P;	    /* ------------ First pass, to get not-precise-enough values for e ------------ */		/* Initialize lengthscale, energy values, and in/decrements at zero */		estimated.l = 0.0;		plus.l      = 0.0;		minus.l     = 0.0;		estimated.thrmv[2] = 0.0;		plus.thrmv[2]      = 0.0;		minus.thrmv[2]     = 0.0;		/* Initialize temperature-dependent values at initial lengthscale */		K0 = K0c/(1.0-estimated.tempv);/*		printf("%17.15g\n",K0); */		Kr = K0;		estimated.thrmv[0] = K0;		plus.thrmv[0]      = K0+DK0;		minus.thrmv[0]     = K0-DK0;		estimated.thrmv[1] = exp(-PISQ*K0*Cc);		plus.thrmv[1]      = exp(-(K0+DK0)*PISQ*Cc);		minus.thrmv[1]     = exp(-(K0-DK0)*PISQ*Cc);/*		printf("%17.15g\t%17.15g\t%17.15g\n",estimated.thrmv[1],plus.thrmv[1],minus.thrmv[1]); */		/* Evolve state-variables in progression to larger length scales via Runge-Kutta method (derk) */		while(estimated.thrmv[0]<Kthresh){			/* Calculate next estimated state and different temperature (K0-value) test states			(to enable a derivative with respect to K0), incrementing lengthscale */			estimated = derk(IN vltRecRel, INOUT estimated, IN dl, IN N);			plus      = derk(IN vltRecRel, INOUT plus,      IN dl, IN N);			minus     = derk(IN vltRecRel, INOUT minus,     IN dl, IN N);			/* Print intermediate numbers to screen (for bug checking)			j++;			if(j%2000==0){				Kr = estimated.thrmv[0]*exp(-estimated.l);				DeDK0plus = (plus.thrmv[2]-estimated.thrmv[2])/DK0;				DeDK0minus = (estimated.thrmv[2]-minus.thrmv[2])/DK0;				DeDK0 = (DeDK0plus+DeDK0minus)/2;				D2eDK02 = (DeDK0plus-DeDK0minus)/DK0;				cp = -R*estimated.thrmv[2]*((6.65e3)/(rho*a0*a0*a0))*T*((T*alpha-2)*alpha-T*DalphaDT+2*(T*alpha-1)*Da0DT/a0+T*(Da0DT*Da0DT)/(a0*a0)-T*D2a0DT2/a0)				  - R*((6.65e3)/(rho*a0*a0*a0))*((T*alpha+2*T*Da0DT/a0)*K0*DeDK0+K0*K0*D2eDK02)				  - T*P*(4.005e-4/rho)*(alpha*alpha-DalphaDT);				printf(         "%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\n",					1-estimated.tempv, estimated.tempv, estimated.l, estimated.thrmv[0], estimated.thrmv[1], estimated.thrmv[2], Kr, Kr/K0, cp);			}			if(j>10000) j=0; */		}		/* Record obtained value of e */		eimprecise = estimated.thrmv[2];	    /* ------------ Second pass, to get precise-enough values for e ------------ */		/* Initialize lengthscale, energy values, and in/decrements at zero */		estimated.l = 0.0;		plus.l      = 0.0;		minus.l     = 0.0;		estimated.thrmv[2] = 0.0;		plus.thrmv[2]      = 0.0;		minus.thrmv[2]     = 0.0;		/* Initialize temperature-dependent values at initial lengthscale */		K0 = K0c/(1.0-estimated.tempv);		Kr = K0;		estimated.thrmv[0] = K0;		plus.thrmv[0]      = K0+DK0;		minus.thrmv[0]     = K0-DK0;		estimated.thrmv[1] = exp(-PISQ*K0*Cc);		plus.thrmv[1]      = exp(-(K0+DK0)*PISQ*Cc);		minus.thrmv[1]     = exp(-(K0-DK0)*PISQ*Cc);		/* Evolve state-variables in progression to larger length scales via Runge-Kutta method (derk) until e stabilizes */		while(estimated.thrmv[0]<Kthresh){			/* Calculate next estimated state and different temperature (K0-value) test states (to enable a derivative with respect to K0), incrementing lengthscale */			estimated = derk(IN vltRecRel, INOUT estimated, IN dl, IN N);			plus      = derk(IN vltRecRel, INOUT plus,      IN dl, IN N);			minus     = derk(IN vltRecRel, INOUT minus,     IN dl, IN N);			/* Check if e has stabilized to some precision (given by digits) and subtract out the first some-odd digits of e by subtracting out eimprecise, and continue working with a relative quantity thrmv[2] = de, where e = eimprecise+de */			if(fabs(estimated.thrmv[2]-eimprecise)<pow(10,-digits)){				estimated.thrmv[2] = estimated.thrmv[2]-eimprecise;				plus.thrmv[2]      = plus.thrmv[2]-eimprecise;				minus.thrmv[2]     = minus.thrmv[2]-eimprecise;				break;			}		}/**/		printf("%17.15g\t%17.15g\t%17.15g\n",estimated.thrmv[2],plus.thrmv[2],minus.thrmv[2]);		/* Continue evolving state-variables in progression to larger length scales via Runge-Kutta method (derk) until K surpasses threshold value Kthresh */		while(estimated.thrmv[0]<Kthresh){			/* Calculate next estimated state and different temperature (K0-value) test states (to enable a derivative with respect to K0), incrementing lengthscale */			estimated = derk(IN vltRecRel, INOUT estimated, IN dl, IN N);			plus      = derk(IN vltRecRel, INOUT plus,      IN dl, IN N);			minus     = derk(IN vltRecRel, INOUT minus,     IN dl, IN N);		}	    /* ------------ By now, should have obtained precise values for e ------------ */		/* Update/set values that depend on evolved state-variables, if not using "Print out intermediate numbers" block of code */		Kr = estimated.thrmv[0]*exp(-estimated.l);		DeDK0plus = (plus.thrmv[2]-estimated.thrmv[2])/DK0;		DeDK0minus = (estimated.thrmv[2]-minus.thrmv[2])/DK0;		DeDK0 = (DeDK0plus+DeDK0minus)/2;		D2eDK02 = (DeDK0plus-DeDK0minus)/DK0;		cp = -R*estimated.thrmv[2]*((6.65e3)/(rho*a0*a0*a0))*T*( (T*alpha-2)*alpha - T*DalphaDT+2*(T*alpha-1)*Da0DT/a0 + T*(Da0DT*Da0DT)/(a0*a0) - T*D2a0DT2/a0 )		   - R*((6.65e3)/(rho*a0*a0*a0))*((T*alpha+2*T*Da0DT/a0)*K0*DeDK0+K0*K0*D2eDK02)		   - T*P*(4.005e-4/rho)*(alpha*alpha-DalphaDT);		/* Print data to screen and output file */		fprintf(outfile,"%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\n",			1-estimated.tempv, estimated.tempv, estimated.l, estimated.thrmv[0], estimated.thrmv[1], eimprecise+estimated.thrmv[2], Kr, Kr/K0, cp);/*		printf(         "%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\n",			1-estimated.tempv, estimated.tempv, estimated.l, estimated.thrmv[0], estimated.thrmv[1], eimprecise+estimated.thrmv[2], Kr, Kr/K0, cp); */	}	/* Close output file */	fclose(outfile);}// vortex loop theory recursion relationsRETARRAY vltRecRel(double l, double z[]){	RETARRAY dzdl;	dzdl.arr[0] = z[0]-A*z[1]*z[0]*z[0];	dzdl.arr[1] = z[1]*(6.0-PISQ*z[0]*(1.0-0.6*log(z[0])));	dzdl.arr[2] = -PI*z[1]*exp(-3.0*l);	return dzdl;}/*  Program Notes  *//*=========================================================================*//*In this program, thrmv[2] sometimes represents e, but sometimes it represents de, where e = eimprecise+de.See further documentation (UCLAresearchNotes.pdf, UCLAresearchPrograms.pdf, derk1.c) for elaboration.plus and minus are test states (with test state-variables)*/