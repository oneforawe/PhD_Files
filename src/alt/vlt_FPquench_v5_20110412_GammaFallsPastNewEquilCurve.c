//  File Comments  //
//=========================================================================//

// I simply deleted references to Kthresh in this version (as well as putting the if-else statement in the FP eqn).  Later, if this version works nicely, I could revert to the previous form of code that doesn't use newlpts.

/* FILENAME: vlt_FPquench.c */
/* VERSION: 5 (2011 Apr 05 - 2011 Apr 12)
            C is now out of the recursion relations, updated fp constants, normal settings */
/* AUTHOR: Andrew Forrester <aforrester@ucla.edu> */
/* DESCRIPTION:
   * 

   * In the basename "vlt_FPquench", "vlt" refers to Vortex Loop Theory (of the superfluid phase transition), "FP" refers to the Fokker-Planck equation that is used in the program, and "quench" refers to the rapid decrease in temperature that occurs in the program calculation.

   Inputs:
   * 

   Output:
   * 
*/
/* EXT FILES: none */
/* COMPILE NOTES:
   * To run, type "g++ -lm vlt_FPquench.c" without the quotes.
   * NOTE: Update P, Cc, K0c, and the output filename before running.
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
const double B     = 4.0*PICU/3.0;
const double THETA = 0.6;

// Program parameters/inputs and data type definitions //
const double lmax     = 10.0;      //  1.0, 10.0  //  maximum of the length scale l (potentially approximate)
const int    lsteps   = 1250;      //  5000, 100000 //  lsteps = lmax/dl = 10/0.0001 = 100,000  (make sure this is divisible by ldatamax)
  //  IF I INCREASE lsteps, I CREATE ERRORS?! (nan)
  //  IT SHOULD GET MORE ACCURATE, NOT WORSE!
const int    lpts     = lsteps+1;  //  from l=0 to l=lmax, inclusive (incl. the boundary-condition-enforcing pnt, there are lpts+1 "l" pnts)
const int    ldatamax = 250;       //  50,100 //  max number of data points recorded per time step
//const double iTfrac   = 0.75;      //  initial temperature fraction Ti/Tkt
//const double qTfrac   = 0.10;      //  quench temperature fraction Tq/Tkt
const double tmax     = 100.0;//10000.0;   //  20600 //  max unitless time (see if(t==...) below)
const double dt0      = 1.0e-4;    //  the time increment (in units of the "diffusion time")
//const double dtf      = 1.0e-4;    //  the time increment (in units of the "diffusion time")
const double a0       = 1.0;       //  a0 in units of a0 (!)
const double a03      = a0*a0*a0;  //  a0 to the third power
const double a06      = a03*a03;   //  a0 to the sixth power
const double FPconst  = 1.0;
const double Lconst   = 1.0;

const double R        = 8.31447215; // Gas constant
const double DK0      = 5e-7;
//const double Kthresh  = 5.0;

// Constants for fpGRecRel //
// taking the derivative of q: "tdq" => "3"
const double fp1 = -3.0+THETA/(1.0-THETA);                        // = -1.50
const double fp2 = (THETA/(1.0-THETA))*(-2.0+THETA/(1.0-THETA));  // = -0.75  not  -3*THETA/(1-THETA) = -4.50
const double fp3 = -2.0+2.0*THETA/(1.0-THETA);                    // =  1.00  not  -3+THETA/(1-THETA) = -1.50
const double fp4 = THETA/(1.0-THETA);                             // =  1.50
/* pulling q out of the derivative: "pqo" => "4"
const double fp1 = -3.0+THETA/(1.0-THETA);  // = -1.50
const double fp2 = -3*THETA/(1-THETA);      // = -4.50
const double fp3 = -3+THETA/(1-THETA);      // = -1.50
const double fp4 = THETA/(1.0-THETA);       // =  1.50
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
double K[lpts+1], G[lpts+1];
double dGdt[1];



//  Function Prototypes  //
//=========================================================================//

// rk4((*f)(),*x,y[],h,n);
// EqRecRel(x,z[],dzdx[],n);
// fpKRecRel(x,Kcalc[1],dKdx[1],n);
// fpGRecRel(x,Dl,Gcalc[lpts+1]);
void rk4(void (*f)(double, double*, double*, unsigned int), double *x, double y[], double h, unsigned n);
void EqRecRel(double x, double z[], double dzdx[], unsigned n);
void fpKRecRel(double x, double Kcalc[1], double dKdx[1], unsigned n);
void fpGRecRel(double x, double Dl, double Gcalc[lpts+1]);



//  Function Definitions  //
//=========================================================================//

main(){
	// Main function definitions //
	int s1,s2;
	double P,Cc,K0c,iTfrac,qTfrac;
	int64_t n_t;
	int j,k;
	double dblsteps=lsteps, dl=lmax/dblsteps;
	double x, l[lpts];
	int    newlpts;
	double K0, y0, z[2];
	double Kr[lpts];
	double Kcalc[lpts+1], Gcalc[lpts+1];
	double k1[lpts], k2[lpts], k3[lpts];
	double testG1[lpts+1], testG2[lpts+1];
	double t, dt, ttrig;//, Dt=dtf-dt0;
	double loopintegral, loopdens;
	FILE *outfile1;  // recording functions versus length-scale l ("vsl"), at selected times
	FILE *outfile2;  // recording functions versus time t ("vst"), at largest length scale lmax
	FILE *outfile3a, *outfile3b;  // recording last set of data to allow resumption of calculations
	char *filename1, *filename2, *filename3a, *filename3b;

	double Tc,T,rho,DrhoDT,D2rhoDT2,Da0DT,D2a0DT2;
	double DeDK0plus,DeDK0minus,DeDK0,D2eDK02,alpha,DalphaDT,cp;
	//double Delok,a;
	// STATE estimated, plus, minus;

	// Pick one pressure (SVP). //Progress through 7 pressures (and associated values of Cc,K0c), where a data file will be produced for each pressure //
	//for(s1=1; s1<=7; s1++){
	//	switch(s1){ // pressure P in bars (See RhosAmps.ods for Cc(P), derived from Ahlers' 1973 data)
			// using normal settings (in vlt_K0cFind.c, vlt_ThermStates.c, and fits) //
			// using normal settings (in vlt_K0cFind.c, vlt_ThermStates.c, and fits) //
	//		case 1:
			P =  0.050;	Cc = 1.10500797392886;	K0c = 0.295614012972922;	//break;
	//		case 2:	P =  1.646;	Cc = 1.0581889918361;	K0c = 0.303871530617906;	break;
	//		case 3:	P =  7.328;	Cc = 0.905918431108949;	K0c = 0.335619263473799;	break;
	//		case 4:	P = 15.031;	Cc = 0.730685115724445;	K0c = 0.385639221938482;	break;
	//		case 5:	P = 18.180;	Cc = 0.667831302441093;	K0c = 0.408941994995917;	break;
	//		case 6:	P = 22.533;	Cc = 0.588281338805611;	K0c = 0.44447352661157;		break;
	//		case 7:	P = 25.868;	Cc = 0.532635362503492;	K0c = 0.474695910004971;	break;
			/* using B' = 2B (in vlt_K0cFind.c, vlt_ThermStates.c, and fits) //
			case 1:	P =  0.050;	Cc = 1.1003376541326;	K0c = 0.338730088317275;
			case 2:	P =  1.646;	Cc = 1.05648724043934;	K0c = 0.348372781936236;	break;
			case 3:	P =  7.328;	Cc = 0.913079877533694;	K0c = 0.385556454984416;	break;
			case 4:	P = 15.031;	Cc = 0.746257959253795;	K0c = 0.444486928182064;	break;
			case 5:	P = 18.180;	Cc = 0.685924754942954;	K0c = 0.472045003545693;	break;
			case 6:	P = 22.533;	Cc = 0.609090754516276;	K0c = 0.51419404658953;		break;
			case 7:	P = 25.868;	Cc = 0.554925448693439;	K0c = 0.550204504048639;	break; */
	//	}

		for(s2=1; s2<=1; s2++){
			switch(s2){
				case 1: iTfrac = 0.999; qTfrac = 0.99; break;
				/*case 2: iTfrac = 0.98; qTfrac = 0.1; break;
				case 3: iTfrac = 0.95; qTfrac = 0.1; break;
				case 4: iTfrac = 0.75; qTfrac = 0.1; break;
				case 5: iTfrac = 0.50; qTfrac = 0.1; break;*/
			}

			// Prepare output files, print identification and values //
			asprintf(&filename1, "vlt_FPquench_vsl_3_T_%g_%g_P_%06.3f_lmax%g_dl%g_tmax%g_dt0%g.dat", iTfrac, qTfrac, P, lmax, dl, tmax, dt0);
			asprintf(&filename2, "vlt_FPquench_vst_3_T_%g_%g_P_%06.3f_lmax%g_dl%g_tmax%g_dt0%g.dat", iTfrac, qTfrac, P, lmax, dl, tmax, dt0);
//			asprintf(&filename3a,"vlt_FPquench_resumeA_T_%g_%g_P_%06.3f_lmax%g_dl%g_tmax%g_dt0%g.dat", iTfrac, qTfrac, P, lmax, dl, tmax, dt0);
//			asprintf(&filename3b,"vlt_FPquench_resumeB_T_%g_%g_P_%06.3f_lmax%g_dl%g_tmax%g_dt0%g.dat", iTfrac, qTfrac, P, lmax, dl, tmax, dt0);
//			outfile3a = fopen(filename3a,"rw");  // E.g., "vlt_FPquench_resumeA.out"
//			outfile3b = fopen(filename3b,"rw");  // E.g., "vlt_FPquench_resumeB.out"
			outfile1  = fopen(filename1,"w");    //  E.g., "vlt_FPquench_vsl.out"
			outfile2  = fopen(filename2,"w");    //  E.g., "vlt_FPquench_vst.out"
			fprintf(outfile1,"# Filename: %s\n", filename1);
			fprintf(outfile2,"# Filename: %s\n", filename2);
			fprintf(outfile1,"# Source: vlt_FPquench.c\n");
			fprintf(outfile2,"# Source: vlt_FPquench.c\n");
			fprintf(outfile1,"# Source version: %s\n", "5 (2011 Apr 05 - ...)");
			fprintf(outfile2,"# Source version: %s\n", "5 (2011 Apr 05 - ...)");
			fprintf(outfile1,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g, B=%21.21g, THETA=%g\n", PI,PISQ,PICU,B,THETA);
			fprintf(outfile2,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g, B=%21.21g, THETA=%g\n", PI,PISQ,PICU,B,THETA);
			fprintf(outfile1,"# Parameter values: iTfrac=%g, qTfrac=%g, P=%g (Cc=%g, K0c=%g), lmax=%g, lsteps=%i, (dl=%g), lpts=%i, ldatamax=%i, tmax=%g, dt0=%g, a0=%g, a03=%g, a06=%g, FPconst=%g, Lconst=%g\n", iTfrac,qTfrac,P,Cc,K0c,lmax,lsteps,dl,lpts,ldatamax,tmax,dt0,a0,a03,a06,THETA,FPconst,Lconst);
			fprintf(outfile2,"# Parameter values: iTfrac=%g, qTfrac=%g, P=%g (Cc=%g, K0c=%g), lmax=%g, lsteps=%i, (dl=%g), lpts=%i, ldatamax=%i, tmax=%g, dt0=%g, a0=%g, a03=%g, a06=%g, FPconst=%g, Lconst=%g\n", iTfrac,qTfrac,P,Cc,K0c,lmax,lsteps,dl,lpts,ldatamax,tmax,dt0,a0,a03,a06,THETA,FPconst,Lconst);

			// Boundary condition enforcement //
			//K[lpts] = G[lpts] = 0;  // K=G=0 at l=lmax+dl, terminates recursion relations at l=lmax
			// instead, use K[lpts] = K[lpts-1], etc. below


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
			x = 0.0;
			l[0] = 0.0;
			K0 = K0c/iTfrac;
			y0 = exp(-PISQ*K0*Cc);
			z[0] = K0;
			z[1] = y0;
			K[0] = z[0];
			Kr[0] = K[0]*exp(-l[0]);
			G[0] = z[1]*exp(-6*l[0])/a06;  // G (not = 0.0)
			loopintegral = 0.0;
			loopintegral += G[0]*exp(3*l[0])*dl;

			// Print initial data to screen and output files //
			printf(          "\n# time step n_t = %ld\tt = %e\t(dt = %e)\n", n_t,t,dt);
			fprintf(outfile1,"\n# time step n_t = %ld\tt = %e\t(dt = %e)\n", n_t,t,dt);
			printf(          "%s\t%s\t%s\t%s\t%s\t%s\t%s\n",   "T/Tc","1-T/Tc","l","K","K/K0","Kr/K0","G");
			fprintf(outfile1,"# %s\t%s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","Kr/K0","G");
			fprintf(outfile2,"# %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "time-step","t","T/Tc","1-T/Tc","l","K","K/K0","Kr/K0","G","loopdens");
			printf(          "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[0], K[0], K[0]/K0, Kr[0]/K0, G[0]);
			fprintf(outfile1,"%g\t%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[0], K[0], K[0]/K0, Kr[0]/K0, G[0]);

			// Calculate state-variables at progressively larger length scales, up to coherence length, via Runge-Kutta method (rk4) //
			i = 0;
			while(i<lpts){
				// Step out in length scale l, calculating next (ith, above zeroth) K, y, l, and G //
				i++;
				rk4(EqRecRel, &x, z, dl, 2);  // equil: rk4 2D (K and y), EqRecRel
				l[i] = l[i-1] + dl;
				K[i] = z[0];
				Kr[i] = K[i]*exp(-l[i]);
				G[i] = z[1]*exp(-6*l[i])/a06;
				loopintegral += G[i]*exp(3*l[i])*dl;
				// Print at most [ldatamax] (e.g., 50) data points to screen and output file 1 //
				if(i%(lsteps/ldatamax)==0){
					printf(          "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, Kr[i]/K0, G[i]);
					fprintf(outfile1,"%g\t%g\t%g\t%g\t%g\t%g\t%g\n", iTfrac, 1-iTfrac, l[i], K[i], K[i]/K0, Kr[i]/K0, G[i]);
				}
			}
			newlpts = i;
			// temporary check:
//			printf("\n\n\n \t\t newlpts = %i\n\n\n\n", newlpts);
			loopdens = 4*PI*a03*loopintegral;
			fflush(outfile1);
			// Boundary, dummy points
			K[newlpts] = K[newlpts-1];
			G[newlpts] = G[newlpts-1];

			// Print out the data at largest length scale to file 2 //
			fprintf(outfile2,"%ld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", n_t, t, iTfrac, 1-iTfrac, l[newlpts-1], K[newlpts-1], K[newlpts-1]/K0, Kr[newlpts-1]/K0, G[newlpts-1], loopdens);
			fflush(outfile2);


			////////////////////////////////////////////////////////////////
			//                                                            //
			//  (n_t>=1) After Instantaneous "Quench" (Temperature Drop)  //
			//                                                            //
			////////////////////////////////////////////////////////////////

			// Initialize //
			K0 = K0c/qTfrac;
			//y0 = exp(-PISQ*K0/2.0);
			K[0] = K0;
			//G[0] = exp(-4.0*l[0])*y0*y0/a04;  // G (not = 0.0)  //  (correct location?)

			// Calculate K at progressively larger length scales, up to coherence length, via Runge-Kutta method (rk4 K) //
			x = 0.0;
			Kcalc[0] = K[0];
			i = 0;
			while(i<newlpts){
				i++;
				rk4(fpKRecRel, &x, Kcalc, dl, 1);  // nequil: rk4 1D (K), fpKRecRel
				K[i] = Kcalc[0];
			}
			newlpts = i;
			// temporary check:
//			printf("\n\n\n \t\t newlpts = %i\n\n\n\n", newlpts);

			// Advance in time: rk3 1D (G), fpGRecRel //
			j = 1;        // for selecting times to print out l-dependent data
			k = 1;        // for selecting times to print out t-dependent data (at largest length scale)
			ttrig = 0.0;  // for selecting
			while(t<=tmax+dt0){  //1e-8){
				n_t++;
//				dt = dt0 + Dt*(1-exp(-t));  // added
				t += dt;
				loopintegral = 0.0;
				printf(          "# time step n_t = %ld\tt = %e\t(dt = %e)\t", n_t,t,dt);
				//fprintf(outfile1,"# time step n_t = %ld\tt = %e\t(dt = %e)\t", n_t,t,dt);

				// Calculate dGdt for all length scales, and update G //
				// nequil: rk4 1D (K), fpKRecRel //
				for(i=1; i<newlpts-1; i++){
					fpGRecRel(l[i], dl, G);
					k1[i]     = dGdt[0]*dt;
					testG1[i] = G[i] + k1[i]/2.0;
				}
				testG1[0]      = testG1[1];
				testG1[newlpts-1] = testG1[newlpts-2];

				for(i=1; i<newlpts-1; i++){
					fpGRecRel(l[i], dl, testG1);
					k2[i]     = dGdt[0]*dt;
					testG2[i] = G[i] - k1[i] + k2[i];
				}
				testG2[0]      = testG2[1];
				testG2[newlpts-1] = testG2[newlpts-2];

				for(i=1; i<newlpts-1; i++){
					fpGRecRel(l[i], dl, testG2);
					k3[i] = dGdt[0]*dt;
					G[i]  = G[i] + k1[i]/6.0 + 2.0*k2[i]/3.0 + k3[i]/6.0;  //  update G here and below**
					loopintegral += G[i]*exp(3*l[i])*dl;
				}
				G[0]      = G[1];       //  **here
				G[newlpts-1] = G[newlpts-2];  //  **and here
				loopintegral += G[0]*exp(3*l[0])*dl;
				loopintegral += G[newlpts-1]*exp(3*l[newlpts-1])*dl;
				loopdens = 4*PI*a03*loopintegral;
				// Boundary, dummy point
				G[newlpts] = G[newlpts-1];

				// Calculate K at progressively larger length scales, up to coherence length, via Runge-Kutta method (rk K) //
				x = 0.0;
				Kcalc[0] = K[0];
				i = 0;
				while(i<newlpts){
					i++;
					rk4(fpKRecRel, &x, Kcalc, dl, 1);  // nequil: rk4 1D (K), fpKRecRel
					K[i] = Kcalc[0];
					Kr[i] = K[i]*exp(-l[i]);
				}
				newlpts = i;
				// temporary check:
//				printf("\n\n\n \t\t newlpts = %i\n\n\n\n", newlpts);
				// Boundary, dummy point
				K[newlpts] = K[newlpts-1];

				// At certain times, print out the results to file 1 //
				//if( (t>=1e-2&&j==1) || (t>=1e-1&&j==2) || (t>=1e0&&j==3) || (t>=2e0&&j==4) || (t>=5e0&&j==5) || (t>=1e1&&j==6) || (t>=2e1&&j==7) || (t>=3e1&&j==8) || (t>=4e1&&j==9) || (t>=5e1&&j==10) ){
				if( (t>=1e-2&&j==1) || (t>=1e-1&&j==2) || (t>=1e0&&j==3) || (t>=1e1&&j==4) || (t>=2e1&&j==5) || (t>=3e1&&j==6) || (t>=4e1&&j==7) || (t>=5e1&&j==8) || (t>=1e2&&j==9) ){  // 1==0    ){// 
					fprintf(outfile1,"\n# time step n_t = %ld\tt = %e\t(dt = %e)\n", n_t,t,dt);
					fprintf(outfile1,"# %s\t%s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","Kr/K0","G");
					for(i=0; i<newlpts; i++){
						// Print at most [ldatamax] (e.g., 1000) data points to screen and output file //
						if(i%(lsteps/ldatamax)==0){
							printf(          "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", qTfrac, 1-qTfrac, l[i], K[i], K[i]/K0, Kr[i]/K0, G[i]);
							fprintf(outfile1,"%g\t%g\t%g\t%g\t%g\t%g\t%g\n", qTfrac, 1-qTfrac, l[i], K[i], K[i]/K0, Kr[i]/K0, G[i]);
						}
					}
					fflush(outfile1);
					j++;
				}

				// Print out the data at the largest length-scale to file 2, selecting data at exponentially sparser time intervals //
				if(t>ttrig){
					fprintf(outfile2,"%ld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", n_t, t, qTfrac, 1-qTfrac, l[newlpts-1], K[newlpts-1], K[newlpts-1]/K0, Kr[newlpts-1]/K0, G[newlpts-1], loopdens);
					fflush(outfile2);
					ttrig += dt0*pow(10.0,0.05*k);
					//ttrig += dt0*exp(0.05*k);
					k++;
				}

				// Bug detection //
				printf(          "%g\t%g\t%g\n", l[newlpts-1], Kr[newlpts-1]/K0, G[newlpts-1]);
				//fprintf(outfile1,"%g\t%g\t%g\n", l[newlpts-1], KG[newlpts-1].c[0]/K0, KG[newlpts-1].c[1]);
			}

			fclose(outfile1);
			fclose(outfile2);
		}
//	}
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
		//f(*x+h*(i/2)/2,s,k[i],n);    // Should this be  f((*x)+h*((i+1)/2)/2,s,k[i],n) ???
		f(*x+h*((i+1)/2)/2,s,k[i],n);  // Correct, right?
	}
	for(j=0;j<n;j++)  y[j] += h*(k[0][j]+2*(k[1][j]+k[2][j])+k[3][j])/6;
	*x += h;
	for(i=0;i<4;i++)  free((char *) k[i]);
	free((char *)s);
}



// Equilibrium (K,y,e) recursion relations //
void EqRecRel(double x, double z[], double dzdx[], unsigned n){
	dzdx[0] = z[0]-B*z[0]*z[0]*z[1];
	if( THETA*log(z[0]) > 0.0 ) // if ac' > a, set ac' = a (set THETA*log(K[i]) to zero)
		dzdx[1] = z[1]*( 6.0 - PISQ*z[0] );
	else
		dzdx[1] = z[1]*( 6.0 - PISQ*z[0]*(1.0-THETA*log(z[0])) );
	//dzdx[2] = -PI*z[1]*exp(-3*x);
}



// Non-equilibrium ("Focker-Planck") K recursion relation //
void fpKRecRel(double x, double Kcalc[1], double dKdx[1], unsigned n){
	dKdx[0] = Kcalc[0]-B*a06*exp(6*x)*Kcalc[0]*Kcalc[0]*G[i-1];
}



// Non-equilibrium Fokker-Planck G recursion relation //
void fpGRecRel(double x, double Dl, double Gcalc[lpts+1]){
	if( THETA*log(K[i]) > 0.0 ) // if ac' > a, set ac' = a (set THETA*log(K[i]) to zero)
		dGdt[0] = 4*FPconst*exp(fp1*x)*(  ( fp2*Gcalc[i] + fp3*(Gcalc[i+1]-Gcalc[i-1])/(2*Dl) + (Gcalc[i+1]-2*Gcalc[i]+Gcalc[i-1])/(Dl*Dl) )/PISQ  +  ( fp4*Gcalc[i]*K[i] + (K[i+1]*Gcalc[i+1]-K[i-1]*Gcalc[i-1])/(2*Dl) )  -  THETA*Gcalc[i]*(K[i+1]-K[i-1])/(2*Dl)  );
	else
		dGdt[0] = 4*FPconst*exp(fp1*x)*(  ( fp2*Gcalc[i] + fp3*(Gcalc[i+1]-Gcalc[i-1])/(2*Dl) + (Gcalc[i+1]-2*Gcalc[i]+Gcalc[i-1])/(Dl*Dl) )/PISQ  +  (1-THETA*log(K[i]))*( fp4*Gcalc[i]*K[i] + (K[i+1]*Gcalc[i+1]-K[i-1]*Gcalc[i-1])/(2*Dl) )  -  THETA*Gcalc[i]*(K[i+1]-K[i-1])/(2*Dl)  );
}



//  Program Notes  //
//=========================================================================//
/*

In this program, thrmv[2] sometimes represents e, but sometimes it represents de, where e = eimprecise+de.

See further documentation (UCLAresearchNotes.pdf, UCLAresearchPrograms.pdf, derk1.c) for elaboration.


plus and minus are test states (with test state-variables)


vortex loop density
loopdens = int_a0^amax G*d^3a = int_a0^amax int_0^2PI G*a^2*dOmega*da = 4PI int_a0^amax G*a^2*da = 
a = a0*exp(l)
da = a0*exp(l)*dl
a^2*da = a0^3*exp(3l)*dl

loopdens = 4PI*a0^3 int_0^lmax G*exp(3l)*dl
loopintegral += G*exp(3l)*dl
loopdens = 4PI*a0^3*loopintegral


vorticity line density
(incorporate fractal dimension) total line-length/volume

linedens = int_a0^amax Gamma*L*
lineintegral


*/
