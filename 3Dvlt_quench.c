//  File Comments  //
//=========================================================================//

/* FILENAME: 3Dvlt_quench.c */
/* VERSION: 10 (2011 Aug 07 - ...)
            With new recursion relations (new B and new y definition with 9 instead of 6) and new THETA.
            With calculated total vortex line density.
            Everywhere with (1-THETA*log(K)) set log(K) to zero if it goes positive.
            With dt determined by limiting relationship between dl and dt.
            With options for 3rd order Runge-Kutta method.
            Got rid of extra point (lpts+1).
            With new large-scale boundary conditions based on number-density-flow derivations.
            With multiple options for boundary conditions of Gamma.
            Adjusted loopintegral by adding 0.5 factor to first point.
            Got rid of temporary dG calculation.
            Bringing back the K,y recursion relations (since the K,G recursion relations produced a bad equilibrium Gamma distribution at Tc -- need more lsteps to make that technique work nicely)
            Not allowing nan's (infinities) at the boundary to propagate (domino-like) over time and blow up Gamma
            Not allowing ac' to exceed a, not allowing Gamma to fall below its equilibrium curve.
            Kthresh is not used since ac' is not allowed to exceed a.
            Now, this version works nicely, so I reverted to the previous form of code that doesn't use newlpts. */
/* AUTHOR: Andrew Forrester <aforrester@ucla.edu> */
/* DESCRIPTION:
   * 

   * In the basename "3Dvlt_quench", "vlt" refers to Vortex Loop Theory (of the superfluid phase transition), "FP" refers to the Fokker-Planck equation that is used in the program, and "quench" refers to the rapid decrease in temperature that occurs in the program calculation.

   Inputs:
   * ...
   * triplets: P,Cc,K0c      (sets of values in 1st "switch" block)
   * couples: iTfrac,qTfrac  (sets of values in 2nd "switch" block)

   Output:
   * 
*/
/* EXT FILES: none */
/* COMPILE NOTES:
   * To compile, type "g++ -lm 3Dvlt_quench.c" without the quotes; then to run, type "./a.out".
   * NOTE: Should try to maximize dt (via dtfrac) before letting program run for a long time.
*/

/* TROUBLESHOOTING IDEAS:
   * Take smaller step sizes (lsteps and time steps dt0)
   * Maybe change Kthresh to something greater than 5.0 ?
*/




//  Function Preparation  //
//=========================================================================//

// Standard routine header files //
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

// Constants //
const double PI    = 3.14159265358979323846;
const double PISQ  = 9.86960440108935799230;
const double PICU  = 31.00627668029981620634;
const double B     = 2.0*PICU*PISQ/3.0;
const double THETA = 0.6;

// Program parameters/inputs and data type definitions //
const int    BC        = 1;          //  Option for Boundary Condition: 1=Hanching's, 2=Gary's, 3=Mine
const double lmax      = 10.0;       //  1.0  //  maximum of the length scale l (potentially approximate)
const int    lsteps    = 5000;      // 12500, 5000, 100000 //  lsteps = lmax/dl = 10/0.0001 = 100,000  (make sure this is divisible by ldatamax)
const int    lpts      = lsteps+1;   //  from l=0 to l=lmax, inclusive
const int    ldatamax  = 100;        //  50,100 //  max number of data points recorded per time step
const double t0        = 1.0;        //  the diffusion time ( a03*(4k/PI)^2/lambda0*kB*T = a0^3*kappa^2*rhos0^2/lambda0*kB*T ) in units of the diffusion time (!)
const double tmax      = 100.0;        //  10000.0 //  max unitless time (see if(t==...) below)
//const double dt0       = 1.0e-4;     //  the time increment (in units of the "diffusion time")
//const double dtf       = 1.0e-4;     //  the time increment (in units of the "diffusion time")
const double a0        = 1.0;        //  a0 in units of a0 (!)
const double a02       = a0*a0;      //  a0 to the second power
const double a03       = a0*a0*a0;   //  a0 to the third power
const double a04       = a02*a02;    //  a0 to the fourth power
const double a06       = a03*a03;    //  a0 to the sixth power
const double Lconst    = 1.0;        //  Constant from L = a0*Lconst*(a/a0)^(1/(1-THETA))  (enters Fokker-Planck eqn via Lambda = lambda0*L)
const double NDflowTOL = 1.0e-6;     //  Number-density flow tolerance (program quits with an error message if numden_out-numden_dec > NDflowTOL)

const double R         = 8.31447215; // Gas constant
const double DK0       = 5e-7;
//const double Kthresh   = 5.0;

// Constants for NEqRecRelG //
// (pulling q out of the derivative)
const double fp1 = THETA/(1.0-THETA)-3.0;   // = -1.50
const double fp2 = THETA/(1.0-THETA);       // =  1.50
const double fp3 = -3.0*THETA/(1.0-THETA);  // = -4.50
const double fp4 = THETA/(1.0-THETA)-3.0;   // = -1.50

// Constant for lineintegral //
const double li = 1.0/(1.0-THETA)+3.0;     // =  5.50

// Parameters for Runge-Kutta method, 3rd order (Version 1) RK3v1
const double B21=0.5,
	     B31=-1.0,    B32=1.0,
	      A1=0.0,      A2=0.0,      A3=0.0,
	      C1=1.0/6.0,  C2=2.0/3.0,  C3=1.0/6.0;
const char*  RKv="RK3v1";
/*
// Parameters for Runge-Kutta method, 3rd order (Version 2) RK3v2
const double B21=0.5,
	     B31=-1.0,    B32=2.0,
	      A1=0.0,      A2=0.5,      A3=1.0,
	      C1=1.0/6.0,  C2=2.0/3.0,  C3=1.0/6.0;
const char*  RKv="RK3v2";

// Parameters for Runge-Kutta method, 3rd order (Version 3) RK3v3
const double B21=1.0/3.0,
	     B31=0.0,     B32=2.0/3.0,
	      A1=0.0,      A2=1.0/3.0,  A3=2.0/3.0,
	      C1=1.0/4.0,  C2=0.0,      C3=3.0/4.0;
const char*  RKv="RK3v3";
*/


// Parameters for Tc(P):
#define Tc0    2.173494256
#define Tc1   -0.009824996
#define Tc2   -0.000118194
#define Tc3   -0.000000437
#define Tc4    0.000000007

// Parameters for rho(T,P):
#define rho00    145.145109496329000
#define rho10   -0.097653969305915
#define rho20    0.334163407001684
#define rho30   -0.446930785976304
#define rho40    0.181879478545246
#define rho01    1.744776044955830
#define rho11   -0.091953899317905
#define rho21    0.179844560873926
#define rho31   -0.133606331352667
#define rho41    0.041022551424992
#define rho02   -0.049165537969017
#define rho12    0.007106988980704
#define rho22   -0.008230542254959
#define rho32    0.000609542602247
#define rho42    0.001149167753923
#define rho03    0.001341503764375
#define rho13   -0.000362007479156
#define rho23    0.000358809384119
#define rho33    0.000064818395436
#define rho43   -0.000104112551303
#define rho04   -0.000016990729415
#define rho14    0.000005538203683
#define rho24   -0.000003157734111
#define rho34   -0.000004999673069
#define rho44    0.000003413312235

// Parameters for alpha(T,P):
#define alpha00   -0.00132545289717059
#define alpha10   0.0179528212646871
#define alpha20   -0.077814417132819
#define alpha30   0.148812649901035
#define alpha40   -0.135348183536229
#define alpha50   0.0575865394848149
#define alpha60   -0.00942356361818271
#define alpha01   0.00169244447357293
#define alpha11   -0.021665108348567
#define alpha21   0.0875997904161722
#define alpha31   -0.155075446681196
#define alpha41   0.133241381828243
#define alpha51   -0.0547135426838175
#define alpha61   0.00857815766443886
#define alpha02   -0.000402016588457985
#define alpha12   0.00500646576912193
#define alpha22   -0.0196997735925578
#define alpha32   0.0343311248462808
#define alpha42   -0.0295791060596021
#define alpha52   0.0123277022988963
#define alpha62   -0.00198704719059122
#define alpha03   3.50822131497139e-05
#define alpha13   -0.000417434943625329
#define alpha23   0.00157100359153513
#define alpha33   -0.0026926366274236
#define alpha43   0.00238134017220718
#define alpha53   -0.00104902823471905
#define alpha63   0.000181757286633977
#define alpha04   -1.51074156982117e-06
#define alpha14   1.69620331490882e-05
#define alpha24   -6.10050371264495e-05
#define alpha34   0.000107369015843627
#define alpha44   -0.000106541517748163
#define alpha54   5.44104698460238e-05
#define alpha64   -1.08691954908159e-05
#define alpha05   3.11237240810798e-08
#define alpha15   -3.20093978259846e-07
#define alpha25   1.08841581501832e-06
#define alpha35   -2.12545135194623e-06
#define alpha45   2.63498247456092e-06
#define alpha55   -1.62910668019171e-06
#define alpha65   3.71525294733195e-07
#define alpha06   -2.37222132465804e-10
#define alpha16   2.08583118201288e-09
#define alpha26   -6.5282210389408e-09
#define alpha36   1.6933030216303e-08
#define alpha46   -2.89786898905826e-08
#define alpha56   2.12469804642866e-08
#define alpha66   -5.29362586788296e-09

// Global variables //
int i;
double K[lpts], G[lpts];
double dGdt[1];



//  Function Prototypes  //
//=========================================================================//

// rk4((*f)(),*x,y[],h,n);
// EqRecRel(x,z[],dzdx[],n);
// NEqRecRelK(x,Kcalc[1],dKdx[1],n);
// NEqRecRelG(x,Dl,Gcalc[lpts+1]);
void rk4(void (*f)(double, double*, double*, unsigned int), double *x, double y[], double h, unsigned n);
void EqRecRel(double x, double z[], double dzdx[], unsigned n);
void NEqRecRelK(double x, double Kcalc[1], double dKdx[1], unsigned n);
void NEqRecRelG(double x, double Dl, double Gcalc[lpts]);



//  Function Definitions  //
//=========================================================================//

int main(){
	// Main function definitions //
	int     s1,s2;
	double  P,Cc,K0c, iTfrac,qTfrac;
	int64_t n_t;
	int     j,k;
	double  dblsteps=lsteps, dl=lmax/dblsteps;
	double  x, l[lpts];
	double  K0, y0, G0final;
	double  Kr[lpts];
	double  z[2], Kcalculate[lpts];
	double  k1[lpts], k2[lpts], k3[lpts];
	double  testG1[lpts], testG2[lpts];
	double  t, dt, dtfrac, ttrig;//, Dt=dtf-dt0;
	double  loopintegral, loopdens, loopdens_old;
	double  lineintegral, linedens;
	double  numden_out, numden_dec;
	FILE    *outfile1;  // recording functions versus length-scale l ("vsl"), at selected times
	FILE    *outfile2;  // recording functions versus time t ("vst"), at largest length scale lmax
	FILE    *outfile3a, *outfile3b;  // recording last set of data to allow resumption of calculations
	char    *filename1, *filename2, *filename3a, *filename3b;

	double  Tc,T,rho,DrhoDT,D2rhoDT2,Da0DT,D2a0DT2;
	double  DeDK0plus,DeDK0minus,DeDK0,D2eDK02,alpha,DalphaDT,cp;
	//double Delok,a;
	// STATE estimated, plus, minus;

	// Time-step size maximized using PDE numerical constraint //
	// Simple constraint:
	dtfrac = 1.0;
	dt = dtfrac*0.5*dl*dl;
	// There could be a more complicated constraint, like Hanching's 2D constraint (I'm not sure where this comes from):
	//dt = 0.9*( 2.0*(dl*dl)/sqrt(4.0/PISQ+iTfrac*K0c*dl*iTfrac*K0c*dl) );

	// Pick one pressure (SVP). //Progress through 7 pressures (and associated values of Cc,K0c), where a data file will be produced for each pressure //
	for(s1=1; s1<=1; s1++){
		switch(s1){ // pressure P in bars (See RhosAmps.ods for Cc(P), derived from Ahlers' 1973 data)
			// using normal settings (in vlt_K0cFind.c, vlt_ThermStates.c, and fits) //
			// using normal settings (in vlt_K0cFind.c, vlt_ThermStates.c, and fits) //
			case 1:	P =  0.050;	Cc = 1.10500797392886;	K0c = 0.295614012972922;	break;/*
			case 2:	P =  1.646;	Cc = 1.0581889918361;	K0c = 0.303871530617906;	break;
			case 3:	P =  7.328;	Cc = 0.905918431108949;	K0c = 0.335619263473799;	break;
			case 4:	P = 15.031;	Cc = 0.730685115724445;	K0c = 0.385639221938482;	break;
			case 5:	P = 18.180;	Cc = 0.667831302441093;	K0c = 0.408941994995917;	break;
			case 6:	P = 22.533;	Cc = 0.588281338805611;	K0c = 0.44447352661157;		break;
			case 7:	P = 25.868;	Cc = 0.532635362503492;	K0c = 0.474695910004971;	break;*/
			/* using B' = 2B (in vlt_K0cFind.c, vlt_ThermStates.c, and fits) //
			case 1:	P =  0.050;	Cc = 1.1003376541326;	K0c = 0.338730088317275;
			case 2:	P =  1.646;	Cc = 1.05648724043934;	K0c = 0.348372781936236;	break;
			case 3:	P =  7.328;	Cc = 0.913079877533694;	K0c = 0.385556454984416;	break;
			case 4:	P = 15.031;	Cc = 0.746257959253795;	K0c = 0.444486928182064;	break;
			case 5:	P = 18.180;	Cc = 0.685924754942954;	K0c = 0.472045003545693;	break;
			case 6:	P = 22.533;	Cc = 0.609090754516276;	K0c = 0.51419404658953;		break;
			case 7:	P = 25.868;	Cc = 0.554925448693439;	K0c = 0.550204504048639;	break; */
		}

		for(s2=1; s2<=1; s2++){
			switch(s2){
				case 1: iTfrac = 0.70; qTfrac = 0.7; break;
				//case 2: iTfrac = 0.98; qTfrac = 0.5; break;
				//case 3: iTfrac = 0.95; qTfrac = 0.5; break;
				//case 4: iTfrac = 0.75; qTfrac = 0.5; break;
				//case 5: iTfrac = 0.50; qTfrac = 0.5; break;
			}

			// Prepare output files, print identification and values //
			asprintf(&filename1, "3Dvlt_quench_vsl_T_%g_%g_P_%06.3f_lmax%g_dl%g_tmax%g_dt%g_BC%i.dat", iTfrac, qTfrac, P, lmax, dl, tmax, dt, BC);
			asprintf(&filename2, "3Dvlt_quench_vst_T_%g_%g_P_%06.3f_lmax%g_dl%g_tmax%g_dt%g_BC%i.dat", iTfrac, qTfrac, P, lmax, dl, tmax, dt, BC);
//			asprintf(&filename3a,"vlt_FPquench_resumeA_T_%g_%g_P_%06.3f_lmax%g_dl%g_tmax%g_dt%g.dat", iTfrac, qTfrac, P, lmax, dl, tmax, dt);
//			asprintf(&filename3b,"vlt_FPquench_resumeB_T_%g_%g_P_%06.3f_lmax%g_dl%g_tmax%g_dt%g.dat", iTfrac, qTfrac, P, lmax, dl, tmax, dt);
//			outfile3a = fopen(filename3a,"rw");  // E.g., "3Dvlt_quench_resumeA.out"
//			outfile3b = fopen(filename3b,"rw");  // E.g., "3Dvlt_quench_resumeB.out"
			outfile1  = fopen(filename1,"w");    //  E.g., "3Dvlt_quench_vsl.out"
			outfile2  = fopen(filename2,"w");    //  E.g., "3Dvlt_quench_vst.out"
			fprintf(outfile1,"# Filename: %s\n", filename1);
			fprintf(outfile2,"# Filename: %s\n", filename2);
			fprintf(outfile1,"# Source: 3Dvlt_quench.c\n");
			fprintf(outfile2,"# Source: 3Dvlt_quench.c\n");
			fprintf(outfile1,"# Source version: %s\n", "10 (2011 Aug 07 - ...)");
			fprintf(outfile2,"# Source version: %s\n", "10 (2011 Aug 07 - ...)");
			fprintf(outfile1,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g, B=%21.21g, THETA=%g\n", PI,PISQ,PICU,B,THETA);
			fprintf(outfile2,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g, B=%21.21g, THETA=%g\n", PI,PISQ,PICU,B,THETA);
			fprintf(outfile1,"# Parameter values: iTfrac=%g, qTfrac=%g, P=%g (Cc=%g, K0c=%g), lmax=%g, lsteps=%i, (dl=%g), lpts=%i, ldatamax=%i, t0=%g, tmax=%g, dt=%g, (dtfrac=%g), a0=%g, a02=%g, a03=%g, a04=%g, a06=%g, Lconst=%g, NDflowTOL=%g, fp1=%g, fp2=%g, fp3=%g, fp4=%g, li=%g\n", iTfrac,qTfrac,P,Cc,K0c,lmax,lsteps,dl,lpts,ldatamax,t0,tmax,dt,dtfrac,a0,a02,a03,a04,a06,Lconst,NDflowTOL,fp1,fp2,fp3,fp4,li);
			fprintf(outfile2,"# Parameter values: iTfrac=%g, qTfrac=%g, P=%g (Cc=%g, K0c=%g), lmax=%g, lsteps=%i, (dl=%g), lpts=%i, ldatamax=%i, t0=%g, tmax=%g, dt=%g, (dtfrac=%g), a0=%g, a02=%g, a03=%g, a04=%g, a06=%g, Lconst=%g, NDflowTOL=%g, fp1=%g, fp2=%g, fp3=%g, fp4=%g, li=%g\n", iTfrac,qTfrac,P,Cc,K0c,lmax,lsteps,dl,lpts,ldatamax,t0,tmax,dt,dtfrac,a0,a02,a03,a04,a06,Lconst,NDflowTOL,fp1,fp2,fp3,fp4,li);
			fprintf(outfile1,"# Option values: BC=%i, RKv=%s\n", BC,RKv); // R,DK0,
			fprintf(outfile2,"# Option values: BC=%i, RKv=%s\n", BC,RKv);


			////////////////////////////////////////////////
			//                                            //
			//  (n_t=0) Initial Equilibrium Calculations  //
			//                                            //
			////////////////////////////////////////////////

			// Initialize time quantities //
			n_t = 0;
			t = 0.0;

			// Initialize quantities at smallest length-scale (l=0) //
			x = l[0] = 0.0;
			z[0]  = K[0] = K0 = K0c/iTfrac;
			z[1]  = y0 = exp(-PISQ*K0*Cc);
			Kr[0] = K[0]*exp(-l[0]);
			G[0]  = z[1]*exp(-6.0*l[0])/a06;
			loopintegral = lineintegral = 0.0;
			loopintegral += 0.5*G[0]*exp(3.0*l[0])*dl;
			lineintegral += 0.5*G[0]*exp(li*l[0])*dl;

			// Print initial data to screen and output files //
			printf(          "\n# time step n_t = %ld\tt = %e\t(dt = %e)\n", n_t,t,dt);
			fprintf(outfile1,"\n# time step n_t = %ld\tt = %e\t(dt = %e)\n", n_t,t,dt);
			printf(            "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","Kr/K0","G");
			fprintf(outfile1,"# %s\t%s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","Kr/K0","G");
			printf(          "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1.0-iTfrac, l[0], K[0], K[0]/K0, Kr[0]/K0, G[0]);
			fprintf(outfile1,"%g\t%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1.0-iTfrac, l[0], K[0], K[0]/K0, Kr[0]/K0, G[0]);
			fprintf(outfile2,"# %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "time-step","t","T/Tc","1-T/Tc","l","K","K/K0","Kr/K0","G","loopdens","linedens");

			// Calculate state-variables at progressively larger length scales, up to coherence length, via Runge-Kutta method (rk4) //
			for(i=1; i<lpts; i++){
				// Step out in length scale l, calculating next (ith, above zeroth) K, y, l, and G //
				rk4(EqRecRel, &x, z, dl, 2);  // equil: rk4 2D (K,y) using EqRecRel
				l[i]  = l[i-1] + dl;
				K[i]  = z[0];
				G[i]  = z[1]*exp(-6.0*l[i])/a06;
				Kr[i] = K[i]*exp(-l[i]);
				loopintegral += G[i]*exp(3.0*l[i])*dl;
				lineintegral += G[i]*exp(li*l[0])*dl;
				// Print at most [ldatamax] (e.g., 50) data points to screen and output file 1 //
				if(i%(lsteps/ldatamax)==0){
					printf(          "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1.0-iTfrac, l[i], K[i], K[i]/K0, Kr[i]/K0, G[i]);
					fprintf(outfile1,"%g\t%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1.0-iTfrac, l[i], K[i], K[i]/K0, Kr[i]/K0, G[i]);
				}
			}
			loopdens = 4.0*PI*a03*loopintegral;
			linedens = 4.0*PI*a04*Lconst*lineintegral;
			fflush(outfile1);

			// Print out the data at largest length scale to file 2 //
			fprintf(outfile2,"%ld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", n_t, t, iTfrac, 1.0-iTfrac, l[lpts-1], K[lpts-1], K[lpts-1]/K0, Kr[lpts-1]/K0, G[lpts-1], loopdens, linedens);
			fflush(outfile2);


			////////////////////////////////////////////////////////////////
			//                                                            //
			//  (n_t>=1) After Instantaneous "Quench" (Temperature Drop)  //
			//                                                            //
			////////////////////////////////////////////////////////////////

			// Initialize //
			K[0] = K0 = K0c/qTfrac;
			Kr[0] = K[0]*exp(-l[0]);

			// Boundary condition: don't let G0 go below this equilibrium value //
			G0final = exp(-6.0*l[0]-PISQ*K0*Cc)/a06;  // = exp(-PISQ*K0*Cc)/a06;

			// Advance in time //
			j = 1;        // for selecting times to print out l-dependent data: outfile1
			k = 1;        // for selecting times to print out t-dependent data (at largest length scale): outfile2
			ttrig = 0.0;  // for selecting times to print out t-dependent data (at largest length scale): outfile2
			//while(t<=tmax+dt){
			while(0==1){
			//while(n_t<=5){
				n_t++;
				t += dt;
				loopdens_old = loopdens;
				loopintegral = lineintegral = 0.0;

				// Bug detection //
				//printf("# time step n_t = %ld\tt = %e\t(dt = %e)\t", n_t,t,dt);
				//printf("# %s\t%s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","Kr/K0","G");

				// Calculate K at progressively larger length scales, up to coherence length, via Runge-Kutta method (rk4 K) //
				x = 0.0;
				Kcalculate[0] = K[0];
				for(i=1; i<lpts; i++){
					rk4(NEqRecRelK, &x, Kcalculate, dl, 1);  // nequil: rk4 1D (K), NEqRecRelK
					K[i] = Kcalculate[0];
					Kr[i] = K[i]*exp(-l[i]);
				}

				// Calculate dGdt for all length scales, and update G, via Runge-Kutta method (rk G) //
				// nequil: rk3 1D (G), NEqRecRelG //
				for(i=1; i<lpts-1; i++){
					NEqRecRelG(l[i]+A1*dl, dl, G);
					k1[i]     = dGdt[0]*dt;
					testG1[i] = G[i] + B21*k1[i];
				}
				k1[0]     = k1[1];                             //*could be questioned
				testG1[0] = G[0] + B21*k1[0];                  //*
				k1[lpts-1]     = 2.0*k1[lpts-2] - k1[lpts-3];  //*
				testG1[lpts-1] = G[lpts-1] + B21*k1[lpts-1];   //*

				for(i=1; i<lpts-1; i++){
					NEqRecRelG(l[i]+A2*dl, dl, testG1);
					k2[i]     = dGdt[0]*dt;
					testG2[i] = G[i] + B31*k1[i] + B32*k2[i];
				}
				k2[0]     = k2[1];                                             //*
				testG2[0] = G[0] + B31*k1[0] + B32*k2[0];                      //*
				k2[lpts-1]     = 2.0*k2[lpts-2] - k2[lpts-3];                  //*
				testG2[lpts-1] = G[lpts-1] + B31*k1[lpts-1] + B32*k2[lpts-1];  //*

				for(i=1; i<lpts-1; i++){
					NEqRecRelG(l[i]+A3*dl, dl, testG2);
					k3[i] = dGdt[0]*dt;
					G[i]  = G[i] + C1*k1[i] + C2*k2[i] + C3*k3[i];
					loopintegral += G[i]*exp(3.0*l[i])*dl;
					lineintegral += G[i]*exp(li*l[i])*dl;
				}
				// Large-scale boundary:
				// no number-density-flow beyond lmax is used to derive:
				// G[lpts-1] = G[lpts-2] + G[lpts-2]*(3.0-PISQ*K[lpts-2]*(1.0-THETA*log(K[lpts-2])))*dl;
				// or G[lpts-1] = G[lpts-2]/( 1.0 - (3.0-PISQ*K[lpts-1]*(1.0-THETA*log(K[lpts-1])))*dl );
				// With the following code, G goes negative (see 3Dvlt_quench_vsl_T_1_0.7_P_00.050_lmax10_dl0.001_tmax1_dt5e-07_BC1_GgoesNeg.dat):
				//if( log(K[lpts-1]) > 0.0 ) // if ac' > a, set ac' = a (replace log(K[lpts-1]) with zero)
				//	G[lpts-1] = G[lpts-2] + G[lpts-2]*(3.0-PISQ*K[lpts-2])*dl;
				//else    G[lpts-1] = G[lpts-2] + G[lpts-2]*(3.0-PISQ*K[lpts-2]*(1.0-THETA*log(K[lpts-2])))*dl;
				// With the following code, G stays positive:
				if( log(K[lpts-1]) > 0.0 ) // if ac' > a, set ac' = a (replace log(K[lpts-1]) with zero)
					G[lpts-1] = G[lpts-2]/( 1.0 - (3.0-PISQ*K[lpts-1])*dl );
				else	G[lpts-1] = G[lpts-2]/( 1.0 - (3.0-PISQ*K[lpts-1]*(1.0-THETA*log(K[lpts-1])))*dl );
				loopintegral += G[lpts-1]*exp(3.0*l[lpts-1])*dl;
				lineintegral += G[lpts-1]*exp(li*l[lpts-1])*dl;
				loopdens = 4.0*PI*a03*loopintegral;  // temporary incomplete calculation of loopdens
				// Small-scale boundary:
				switch(BC){
					case 1: // Hanching's bndry cond. Second l-derivative of G is set to zero (linear continuation).
						// best BC choice since it's simple and BC3 is checked to be true below (numden_out=numden_dec)
						G[0] = 2.0*G[1] - G[2];
						if(G[0]<=G0final)  G[0] = G0final;
						break;
					case 2: // Gary's bndry cond. Immediate and permanent drop to new equilibrium value.
						G[0] = G0final;
						break;
					case 3: // My bndry cond. Total number-density outflow must occur at smallest scale.
						if( log(K[0]) > 0.0 ) // if ac' > a, set ac' = a (replace log(K[lpts-1]) with zero)
							G[0] = (PI*t0*(loopdens_old-loopdens)*dl/(16.0*Lconst*a03) - G[1]*dt) / ( -(3.0-PISQ*K[0])*dl*dt + PISQ*t0*dl*dl/(8.0*Lconst) - dt);
						else	G[0] = (PI*t0*(loopdens_old-loopdens)*dl/(16.0*Lconst*a03) - G[1]*dt) / ( -(3.0-PISQ*K[0]*(1.0-THETA*log(K[0])))*dl*dt + PISQ*t0*dl*dl/(8.0*Lconst) - dt);
						if(G[0]<=G0final)  G[0] = G0final;
						break;
				}
				loopintegral += 0.5*G[0]*exp(3.0*l[0])*dl;
				lineintegral += 0.5*G[0]*exp(li*l[0])*dl;
				loopdens = 4.0*PI*a03*loopintegral;
				linedens = 4.0*PI*a04*Lconst*lineintegral;

				// Check to make sure number density outflow matches number density decrease //
				if( log(K[0]) > 0.0 ) // if ac' > a, set ac' = a (replace log(K[lpts-1]) with zero)
					numden_out = 4.0*Lconst*(a0/t0)*( G[0]*K[0] + (-3.0*G[0]+(G[1]-G[0])/dl)/PISQ )*4.0*PI*a02;
				else	numden_out = 4.0*Lconst*(a0/t0)*( G[0]*K[0]*(1.0-THETA*log(K[0])) + (-3.0*G[0]+(G[1]-G[0])/dl)/PISQ )*4.0*PI*a02;
				numden_dec = -(loopdens-loopdens_old)/dt;
				if(numden_out-numden_dec > NDflowTOL){
					printf("\n\nnumden_out - numden_dec = %g", numden_out-numden_dec);
					printf("\nCode Problem: Vortex-loop number-density outflow doesn't equal decrease in number-density!");
					printf("\n(This probably means that the time-step size dt is too large.)\n");
					exit(EXIT_FAILURE);
				}

				// At certain times, print out the results to file 1 //
				//if( (t>=1e-2&&j==1) || (t>=1e-1&&j==2) || (t>=1e0&&j==3) || (t>=1e1&&j==4) || (t>=2e1&&j==5) || (t>=3e1&&j==6) || (t>=4e1&&j==7) || (t>=5e1&&j==8) || (t>=1e2&&j==9) ){
				//if( (t>=1e-4&&j==1) || (t>=1e-3&&j==2) || (t>=1e-2&&j==3) || (t>=1e-1&&j==4) || (t>=1e0&&j==5) || (t>=1e1&&j==6) || (t>=2e1&&j==7) || (t>=3e1&&j==8) || (t>=4e1&&j==9) || (t>=5e1&&j==10) || (t>=1e2&&j==11) ){
				if( (t>=1e-4&&j==1) || (t>=1e-3&&j==2) || (t>=1e-2&&j==3) || (t>=1e-1&&j==4) || (t>=2e-1&&j==5) || (t>=3e-1&&j==6) || (t>=4e-1&&j==7) || (t>=5e-1&&j==8) || (t>=6e-1&&j==9) || (t>=7e-1&&j==10) || (t>=8e-1&&j==11) || (t>=9e-1&&j==12) || (t>=1e0&&j==13) ){
				//if( n_t==1 || n_t==2 || n_t==3 || n_t==4 || n_t==5 || n_t==6 || n_t==7 || n_t==8 || n_t==9 || n_t==10 ){  // 1==0    ){// 
					fprintf(outfile1,"\n# time step n_t = %ld\tt = %e\t(dt = %e)\n", n_t,t,dt);
					fprintf(outfile1,"# %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","Kr/K0","G","dG");
					for(i=0; i<lpts; i++){
						// Print at most [ldatamax] (e.g., 1000) data points to screen and output file //
						if(i%(lsteps/ldatamax)==0){
							printf(          "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", qTfrac, 1.0-qTfrac, l[i], K[i], K[i]/K0, Kr[i]/K0, G[i]);
							fprintf(outfile1,"%g\t%g\t%g\t%g\t%g\t%g\t%g\n", qTfrac, 1.0-qTfrac, l[i], K[i], K[i]/K0, Kr[i]/K0, G[i]);
						}
					}
					fflush(outfile1);
					j++;
				}

				// Print out the data at the largest length-scale to file 2, selecting data at exponentially sparser time intervals //
				if(t>ttrig){
					fprintf(outfile2,"%ld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", n_t, t, qTfrac, 1.0-qTfrac, l[lpts-1], K[lpts-1], K[lpts-1]/K0, Kr[lpts-1]/K0, G[lpts-1], loopdens, linedens);
					fflush(outfile2);
					ttrig += dt*pow(10.0,0.05*k);
					//ttrig += dt0*exp(0.05*k);
					k++;
					// Sparse update //
					printf(          "n_t = %ld  t = %e  (dt = %e)  l = %g  Kr/K0 = %g  G = %g  loopdens = %g  linedens = %g  (Tfrac: %g->%g)\n", n_t, t, dt, l[lpts-1], Kr[lpts-1]/K0, G[lpts-1], loopdens, linedens, iTfrac, qTfrac);
				}

				// Bug detection //
				//printf(          "%g\t%g\t%g\n", l[lpts-1], Kr[lpts-1]/K0, G[lpts-1]);
				//fprintf(outfile1,"%g\t%g\t%g\n", l[lpts-1], KG[lpts-1].c[0]/K0, KG[lpts-1].c[1]);
			}

			fclose(outfile1);
			fclose(outfile2);
		}
	}
}








/*
		// Loop through different temperatures (calculations for each temp are independent) //
		//for(i=0; i<35; i++){
		for(i=0; i<47; i++){
			// Decrement temperature (further below Tc) //
			estimated.tempv = (1.0e-8)*exp(0.4*i);

			// Set temperature- and pressure-dependent values //
			Tc = Tc0 + Tc1*P + Tc2*P*P + Tc3*P*P*P + Tc4*P*P*P*P;
			T = Tc*(1-estimated.tempv);
			rho = rho00         + rho10*T         + rho20*T*T         + rho30*T*T*T         + rho40*T*T*T*T 
			    + rho01*P       + rho11*T*P       + rho21*T*T*P       + rho31*T*T*T*P       + rho41*T*T*T*T*P 
			    + rho02*P*P     + rho12*T*P*P     + rho22*T*T*P*P     + rho32*T*T*T*P*P     + rho42*T*T*T*T*P*P 
			    + rho03*P*P*P   + rho13*T*P*P*P   + rho23*T*T*P*P*P   + rho33*T*T*T*P*P*P   + rho43*T*T*T*T*P*P*P 
			    + rho04*P*P*P*P + rho14*T*P*P*P*P + rho24*T*T*P*P*P*P + rho34*T*T*T*P*P*P*P + rho44*T*T*T*T*P*P*P*P;
			DrhoDT = rho10         + 2*rho20*T         + 3*rho30*T*T         + 4*rho40*T*T*T 
			       + rho11*P       + 2*rho21*T*P       + 3*rho31*T*T*P       + 4*rho41*T*T*T*P 
			       + rho12*P*P     + 2*rho22*T*P*P     + 3*rho32*T*T*P*P     + 4*rho42*T*T*T*P*P 
			       + rho13*P*P*P   + 2*rho23*T*P*P*P   + 3*rho33*T*T*P*P*P   + 4*rho43*T*T*T*P*P*P 
			       + rho14*P*P*P*P + 2*rho24*T*P*P*P*P + 3*rho34*T*T*P*P*P*P + 4*rho44*T*T*T*P*P*P*P;
			D2rhoDT2 = 2*rho20         + 6*rho30*T         + 12*rho40*T*T 
			         + 2*rho21*P       + 6*rho31*T*P       + 12*rho41*T*T*P 
			         + 2*rho22*P*P     + 6*rho32*T*P*P     + 12*rho42*T*T*P*P 
			         + 2*rho23*P*P*P   + 6*rho33*T*P*P*P   + 12*rho43*T*T*P*P*P 
			         + 2*rho24*P*P*P*P + 6*rho34*T*P*P*P*P + 12*rho44*T*T*P*P*P*P;
			a0 = 549*Tc*K0c/rho;  // in angstroms
			// a0(T,P) = ( m^2 kB Tc(P) K0c(P) ) / (hbar^2 rho(T,P) )
			// where m = 6.65e-27 is the mass of 4He atom in kg, kB is the Boltzmann constant, and hbar is the Dirac constant
			Da0DT = -549*Tc*K0c*DrhoDT/(rho*rho);
			D2a0DT2 = 2*549*Tc*K0c*DrhoDT*DrhoDT/(rho*rho*rho) - 549*Tc*K0c*D2rhoDT2/(rho*rho);
			alpha = alpha00             + alpha10*T             + alpha20*T*T             + alpha30*T*T*T             + alpha40*T*T*T*T             + alpha50*T*T*T*T*T             + alpha60*T*T*T*T*T*T 
			      + alpha01*P           + alpha11*T*P           + alpha21*T*T*P           + alpha31*T*T*T*P           + alpha41*T*T*T*T*P           + alpha51*T*T*T*T*T*P           + alpha61*T*T*T*T*T*T*P 
			      + alpha02*P*P         + alpha12*T*P*P         + alpha22*T*T*P*P         + alpha32*T*T*T*P*P         + alpha42*T*T*T*T*P*P         + alpha52*T*T*T*T*T*P*P         + alpha62*T*T*T*T*T*T*P*P 
			      + alpha03*P*P*P       + alpha13*T*P*P*P       + alpha23*T*T*P*P*P       + alpha33*T*T*T*P*P*P       + alpha43*T*T*T*T*P*P*P       + alpha53*T*T*T*T*T*P*P*P       + alpha63*T*T*T*T*T*T*P*P*P 
			      + alpha04*P*P*P*P     + alpha14*T*P*P*P*P     + alpha24*T*T*P*P*P*P     + alpha34*T*T*T*P*P*P*P     + alpha44*T*T*T*T*P*P*P*P     + alpha54*T*T*T*T*T*P*P*P*P     + alpha64*T*T*T*T*T*T*P*P*P*P 
			      + alpha05*P*P*P*P*P   + alpha15*T*P*P*P*P*P   + alpha25*T*T*P*P*P*P*P   + alpha35*T*T*T*P*P*P*P*P   + alpha45*T*T*T*T*P*P*P*P*P   + alpha55*T*T*T*T*T*P*P*P*P*P   + alpha65*T*T*T*T*T*T*P*P*P*P*P 
			      + alpha06*P*P*P*P*P*P + alpha16*T*P*P*P*P*P*P + alpha26*T*T*P*P*P*P*P*P + alpha36*T*T*T*P*P*P*P*P*P + alpha46*T*T*T*T*P*P*P*P*P*P + alpha56*T*T*T*T*T*P*P*P*P*P*P + alpha66*T*T*T*T*T*T*P*P*P*P*P*P;
			DalphaDT = alpha10             + 2*alpha20*T             + 3*alpha30*T*T             + 4*alpha40*T*T*T             + 5*alpha50*T*T*T*T             + 6*alpha60*T*T*T*T*T 
			         + alpha11*P           + 2*alpha21*T*P           + 3*alpha31*T*T*P           + 4*alpha41*T*T*T*P           + 5*alpha51*T*T*T*T*P           + 6*alpha61*T*T*T*T*T*P 
			         + alpha12*P*P         + 2*alpha22*T*P*P         + 3*alpha32*T*T*P*P         + 4*alpha42*T*T*T*P*P         + 5*alpha52*T*T*T*T*P*P         + 6*alpha62*T*T*T*T*T*P*P 
			         + alpha13*P*P*P       + 2*alpha23*T*P*P*P       + 3*alpha33*T*T*P*P*P       + 4*alpha43*T*T*T*P*P*P       + 5*alpha53*T*T*T*T*P*P*P       + 6*alpha63*T*T*T*T*T*P*P*P 
			         + alpha14*P*P*P*P     + 2*alpha24*T*P*P*P*P     + 3*alpha34*T*T*P*P*P*P     + 4*alpha44*T*T*T*P*P*P*P     + 5*alpha54*T*T*T*T*P*P*P*P     + 6*alpha64*T*T*T*T*T*P*P*P*P 
			         + alpha15*P*P*P*P*P   + 2*alpha25*T*P*P*P*P*P   + 3*alpha35*T*T*P*P*P*P*P   + 4*alpha45*T*T*T*P*P*P*P*P   + 5*alpha55*T*T*T*T*P*P*P*P*P   + 6*alpha65*T*T*T*T*T*P*P*P*P*P 
			         + alpha16*P*P*P*P*P*P + 2*alpha26*T*P*P*P*P*P*P + 3*alpha36*T*T*P*P*P*P*P*P + 4*alpha46*T*T*T*P*P*P*P*P*P + 5*alpha56*T*T*T*T*P*P*P*P*P*P + 6*alpha66*T*T*T*T*T*P*P*P*P*P*P;

			// Initialize lengthscale, energy values, and in/decrements at zero //
			estimated.l = 0.0;
			plus.l      = 0.0;
			minus.l     = 0.0;
			estimated.thrmv[2] = 0.0;
			plus.thrmv[2]      = 0.0;
			minus.thrmv[2]     = 0.0;

			// Initialize temperature-dependent values at initial lengthscale //
			K0 = K0c/(1.0-estimated.tempv);
			Kr = K0;
			estimated.thrmv[0] = K0;
			plus.thrmv[0]      = K0+DK0;
			minus.thrmv[0]     = K0-DK0;
			estimated.thrmv[1] = exp(-PISQ*K0*Cc);        //  exp(-PISQ*K0*Cc)        1.0/(exp(PISQ*K0*Cc)-1.0)
			plus.thrmv[1]      = exp(-(K0+DK0)*PISQ*Cc);  //  exp(-(K0+DK0)*PISQ*Cc)  1.0/(exp((K0+DK0)*PISQ*Cc)-1.0)
			minus.thrmv[1]     = exp(-(K0-DK0)*PISQ*Cc);  //  exp(-(K0-DK0)*PISQ*Cc)  1.0/(exp((K0-DK0)*PISQ*Cc)-1.0)

			printf("%17.15g\t%17.15g\t%17.15g\n", estimated.thrmv[1], plus.thrmv[1], minus.thrmv[1]);

			// Evolve state-variables in progression to larger length scales via Runge-Kutta method (rk4) //
			while(estimated.thrmv[0]<Kthresh){
				// Calculate next estimated state and different temperature (K0-value) test states 
				// (to enable a derivative with respect to K0), incrementing lengthscale //
				estimated = derk(IN vltRecRel, INOUT estimated, IN dl, IN N);
				plus      = derk(IN vltRecRel, INOUT plus,      IN dl, IN N);
				minus     = derk(IN vltRecRel, INOUT minus,     IN dl, IN N);
				// Print intermediate numbers to screen (for bug checking)
				//j++;
				//if(j%2000==0){
				//	Kr = estimated.thrmv[0]*exp(-estimated.l);
				//	DeDK0plus = (plus.thrmv[2]-estimated.thrmv[2])/DK0;
				//	DeDK0minus = (estimated.thrmv[2]-minus.thrmv[2])/DK0;
				//	DeDK0 = (DeDK0plus+DeDK0minus)/2;
				//	D2eDK02 = (DeDK0plus-DeDK0minus)/DK0;
				//	cp = -R*estimated.thrmv[2]*((6.65e3)/(rho*a0*a0*a0))*T*((T*alpha-2)*alpha-T*DalphaDT+2*(T*alpha-1)*Da0DT/a0+T*(Da0DT*Da0DT)/(a0*a0)-T*D2a0DT2/a0)
				//	  - R*((6.65e3)/(rho*a0*a0*a0))*((T*alpha+2*T*Da0DT/a0)*K0*DeDK0+K0*K0*D2eDK02)
				//	  - T*P*(4.005e-4/rho)*(alpha*alpha-DalphaDT);
				//	printf(         "%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\n",
				//		1-estimated.tempv, estimated.tempv, estimated.l, estimated.thrmv[0], estimated.thrmv[1], estimated.thrmv[2], Kr, Kr/K0, cp);
				//}
				//if(j>10000) j=0;
			}

			// Update/set values that depend on evolved state-variables, if not using "Print intermediate numbers..." block of code //
			Kr = estimated.thrmv[0]*exp(-estimated.l);
			DeDK0plus = (plus.thrmv[2]-estimated.thrmv[2])/DK0;
			DeDK0minus = (estimated.thrmv[2]-minus.thrmv[2])/DK0;
			DeDK0 = (DeDK0plus+DeDK0minus)/2;
			D2eDK02 = (DeDK0plus-DeDK0minus)/DK0;
			cp = -R*estimated.thrmv[2]*((6.65e3)/(rho*a0*a0*a0))*T*( (T*alpha-2)*alpha - T*DalphaDT + 2*(T*alpha-1)*Da0DT/a0 + T*(Da0DT*Da0DT)/(a0*a0) - T*D2a0DT2/a0 )
			   - R*((6.65e3)/(rho*a0*a0*a0))*((T*alpha+2*T*Da0DT/a0)*K0*DeDK0+K0*K0*D2eDK02)
			   - T*P*(4.005e-4/rho)*(alpha*alpha-DalphaDT);

			// Print data to screen and output file //
			fprintf(outfile,"%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\n",
				1-estimated.tempv, estimated.tempv, estimated.l, estimated.thrmv[0], estimated.thrmv[1], estimated.thrmv[2], Kr, Kr/K0, cp);
			printf(         "%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\n",
				1-estimated.tempv, estimated.tempv, estimated.l, estimated.thrmv[0], estimated.thrmv[1], estimated.thrmv[2], Kr, Kr/K0, cp);
		}

		// Close output file //
		fclose(outfile);
	}
}
*/



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



// Equilibrium (K,y,e) recursion relations //
void EqRecRel(double x, double z[], double dzdx[], unsigned n){
	dzdx[0] = z[0]-B*z[0]*z[0]*z[1];
	if( log(z[0]) > 0.0 ) // if ac' > a, set ac' = a (replace log(K[i]) with zero)
		dzdx[1] = z[1]*( 9.0 - PISQ*z[0] );
	else
		dzdx[1] = z[1]*( 9.0 - PISQ*z[0]*(1.0-THETA*log(z[0])) );
	//dzdx[2] = -PI*z[1]*exp(-3*x);
}



// Non-equilibrium K recursion relation //
void NEqRecRelK(double x, double Kcalc[1], double dKdx[1], unsigned n){
	dKdx[0] = Kcalc[0]-0.5*B*a06*exp(6.0*x)*G[i-1]*Kcalc[0]*Kcalc[0];
}



// Non-equilibrium Fokker-Planck G recursion relation //
void NEqRecRelG(double x, double Dl, double Gcalc[lpts+1]){
	if( log(K[i]) > 0.0 ) // if ac' > a, set ac' = a (replace log(K[i]) with zero)
		dGdt[0] = 4.0*Lconst*exp(fp1*x)*(  ( fp2*Gcalc[i]*K[i] + (K[i+1]*Gcalc[i+1]-K[i-1]*Gcalc[i-1])/(2.0*Dl) )  -  THETA*Gcalc[i]*(K[i+1]-K[i-1])/(2.0*Dl)  +  ( fp3*Gcalc[i] + fp4*(Gcalc[i+1]-Gcalc[i-1])/(2.0*Dl) + (Gcalc[i+1]-2.0*Gcalc[i]+Gcalc[i-1])/(Dl*Dl) )/PISQ  );
	else
		dGdt[0] = 4.0*Lconst*exp(fp1*x)*(  (1.0-THETA*log(K[i]))*( fp2*Gcalc[i]*K[i] + (K[i+1]*Gcalc[i+1]-K[i-1]*Gcalc[i-1])/(2.0*Dl) )  -  THETA*Gcalc[i]*(K[i+1]-K[i-1])/(2.0*Dl)  +  ( fp3*Gcalc[i] + fp4*(Gcalc[i+1]-Gcalc[i-1])/(2.0*Dl) + (Gcalc[i+1]-2.0*Gcalc[i]+Gcalc[i-1])/(Dl*Dl) )/PISQ  );
}



//  Program Notes  //
//=========================================================================//
/*

See further documentation (UCLAresearchNotes.pdf, UCLAresearchPrograms.pdf, derk1.c) for elaboration.
plus and minus are test states (with test state-variables)



*/
