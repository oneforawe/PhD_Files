//  File Comments  //
//=========================================================================//

/* FILENAME: 3Dvlt_macro_cp.c */
/* VERSION: 9 (2011 Dec 14 - 2012 Jul 31)
            Added temperature Op (Options), similar to other 2Dvpt and 3Dvlt programs.
            Got rid of alpha, since it was not appropriate for this application (i.e., it wasn't experimental data - it was incomplete theoretical data).
            Added more output quantities.
            More terms have been added, due to the refined derivation of cp.
            New name: from vlt_HeatCap_P.c to 3Dvlt_macro_cp.c
            C is now out of the recursion relations, normal settings */
/* AUTHOR: Andrew Forrester <aforrester@ucla.edu> */
/* DESCRIPTION:
   * This C program calculates the (length-scale-dependent) molar specific heat capacity at constant pressure c_P (cp) (and other thermodynamic quantities) of a system of liquid helium-4 (4He) at the macroscopic length scale, at various temperatures and pressures.
   * The system is a spherical region (of variable diameter) of liquid 4He at some constant pressure, near or below the critical transition temperature Tc (which depends on pressure) separating normal fluid and superfluid (aka the lambda temperature).
   * The program uses the vortex-loop theory of the superfluid phase transition (i.e., its differential recursion relations) of Gary A. Williams and Subodh R. Shenoy and a 4th-order Runge-Kutta integration method to calculate (length-scale-dependent) thermodynamic quantities (especially cp) of the system at a set of temperatures {Tn}, at various length scales (no greater than a length scale corresponding to the the limit K<Kthresh), and at some pressure (which also determines the choice of Cc, the critical dimensionless vortex-loop core parameter, and K0c, the critical "bare" coupling).

   * Explain there are 7 pressures that are explored...

   * In the basename "vlt_HeatCap_P", "vlt" refers to Vortex Loop Theory (of the superfluid phase transition), "HeatCap" refers to and the molar specific HEAT CAPacity of the system under examination, and "P" refers to the pressure, which must be specified in units of bars in the parameters/inputs below.
   * An alternate of this program, using mpfr, is vlt_pHeatCap_P.c.

   * This program could be subsumed in vlt_ThermStates.c...

   Inputs:
   * Dexp    - Why not?
   * dl      - step in l, the length scale, which allows derivatives with respect to l
   * DK0     - step in K0, the "bare" dimensionless superfluid ratio, which allows derivatives with respect to K0 (corresponds to derivatives wrt temperature)
   * Kthresh - (... improve this later: K0c is calculated to be critical within a certain length scale region; as soon as K diverges, you are out of that region and should stop assuming that K0c is truely the critical value of K0. Is this right?)
   * You could add options for selecting different sets of temperatures.  Right now, you can just "Decrement temperature (further below Tc)".
   * 7 sets of 3-tuples (P,Cc,K0c) from vlt_K0cFind.c and analysis using Ahlers' 1973 data.
   * Parameters for Tc(P), rho(T,P), and alpha(T,P).

   Output:
   * This program returns text file of data in a two dimensional array.
*/
/* EXT FILES: vlt_derk.c */
/* COMPILE NOTES:
   * To run, type "g++ -lm 3Dvlt_macro_cp.c vlt_derk.c" without the quotes.
   * NOTE: Be sure that N=3 in vlt_derk.c when compiling.
   * NOTE: Update P, Cc, K0c, and the output filename before running.
*/


//  Function Preparation  //
//=========================================================================//

// Standard routine header files //
#include <stdio.h>
#include <math.h>

// Constants and labels //
#define IN           // input label
#define OUT          // output label
#define INOUT        // input/output label
const double PI    = 3.14159265358979323846;
const double PISQ  = 9.86960440108935861883;
const double PICU  = 31.00627668029981620634;
const double B     = 4.0*PICU/3.0;
const double THETA = 0.6;
const double R     = 8.31447215;

// Program parameters/inputs //
//#define Dexp, lmax?
const int    Op      = 7;        //  Option (1 -> hand-picked temp's; else temperature-spread formulas)
const int    N       = 3;
const double Kthresh = 5.0;
const double DK0     = 5e-7;
const double dl      = 0.0001;
const int    lmax    = 10;       // 100 corresponds to largest length scale D = a0*exp(lmax) ~ 6.8e33 meters (huge diameter system)

// Parameters for Tc(P):
const double Tc0 = 2.17349425585161;
const double Tc1 = -0.00982499579394534;
const double Tc2 = -0.000118194448444384;
const double Tc3 = -4.36914591522034e-07;
const double Tc4 = 7.39407378262721e-09;

// Parameters for rho(T,P):
const double rho00 = 145.145109496329;
const double rho10 = -0.0976539693059151;
const double rho20 = 0.334163407001684;
const double rho30 = -0.446930785976304;
const double rho40 = 0.181879478545246;
const double rho01 = 1.74477604495583;
const double rho11 = -0.0919538993179052;
const double rho21 = 0.179844560873926;
const double rho31 = -0.133606331352667;
const double rho41 = 0.0410225514249919;
const double rho02 = -0.0491655379690169;
const double rho12 = 0.00710698898070406;
const double rho22 = -0.00823054225495917;
const double rho32 = 0.000609542602247143;
const double rho42 = 0.00114916775392305;
const double rho03 = 0.0013415037643754;
const double rho13 = -0.000362007479155809;
const double rho23 = 0.000358809384119286;
const double rho33 = 6.48183954357527e-05;
const double rho43 = -0.000104112551302631;
const double rho04 = -1.69907294147191e-05;
const double rho14 = 5.53820368251513e-06;
const double rho24 = -3.15773411117433e-06;
const double rho34 = -4.99967306908062e-06;
const double rho44 = 3.41331223468399e-06;

// Data type definitions //
typedef struct {
	double l;
	double tempv;
	double thrmv[N];
} STATE;

typedef struct {
	double arr[N];
} RETARRAY;



//  Function Prototypes  //
//=========================================================================//

//   derk(func(),calculated,Dl,n)
//   vltRecRel(l,z)
extern STATE derk(RETARRAY (*func)(double l, double z[]), STATE calculated, double Dl, unsigned n);
RETARRAY vltRecRel(double l, double z[]);



//  Function Definitions  //
//=========================================================================//

main(){
	// Main function definitions //
	int    s,i,j=0, notdone=true, TempLabelOp;
	double Tfrac, Tempv;
	double P,Cc,K0c;
	double K0,Kr,Tc,T,rho,DrhoDT,D2rhoDT2,a0,Da0DT,D2a0DT2;
	double DeDK0plus,DeDK0minus,DeDK0,D2eDK02,alpha,DalphaDT,cp;
	double cp01,cp02,cp03,cp04,cp05,cp06,cp07;
	double amax,G;
	STATE estimated, plus, minus;
	FILE *outfile;
	char *filename;

	switch(Op){  //  TempLabelOp -> 0 or 1 for labeling the file with "T" or "Tv", respectively
		case 1:  TempLabelOp=0;  break;
		case 2:  TempLabelOp=1;  break;
		case 3:  TempLabelOp=0;  break;
		case 4:  TempLabelOp=0;  break;
		case 5:  TempLabelOp=1;  break;
		case 6:  TempLabelOp=1;  break;
		case 7:  TempLabelOp=0;  break;
	} // TempLabelOp, and the file label (T or Tv) is meant to indicate which quantity is best used for plotting the data T/Tc = T/Tkt or Tv = tempv = 1-T/Tc = 1-T/Tkt

	// Progress through 7 pressures (and associated values of Cc,K0c), where a data file will be produced for each pressure //
	//for(s=1; s<=7; s++){
	//for(s=1; s<=6; s++){
	for(s=1; s<=5; s++){
		switch(s){ // pressure P in bars
			/* using normal settings (with "noifthen" in 3Dvlt_K0cFind.c v9 Output30, 3Dvlt_macro_noifthen.c,
			// and fits in results/plots_macro/3D_12_Thesis/) //
			case 1:	P =  0.050;	Cc = 1.10500984397513;	K0c = 0.295613694950781;	break;
			case 2:	P =  1.646;	Cc = 1.05819110517992;	K0c = 0.303871144015907;	break;
			case 3:	P =  7.328;	Cc = 0.905919737570672;	K0c = 0.335618952610176;	break;
			case 4:	P = 15.031;	Cc = 0.730687687879284;	K0c = 0.385638338959041;	break;
			case 5:	P = 18.180;	Cc = 0.667832798147861;	K0c = 0.408941395856842;	break;
			case 6:	P = 22.533;	Cc = 0.588282353648267;	K0c = 0.444473020755711;	break;
			case 7:	P = 25.868;	Cc = 0.532636580468056;	K0c = 0.474695188872846;	break;*/

			// using normal settings (with "noifthen" in 3Dvlt_K0cFind.c v9 Output33, 3Dvlt_macro_noifthen.c,
			// and fits in results/plots_macro/3D_12_Thesis/) //
			//case 1:	P =  0.05;	Cc = 1.10500984397513;	K0c = 0.295613694950781;	break;
			//case 2:	P =  7.27;	Cc = 0.907367789194129;	K0c = 0.335274865837589;	break;
			//case 3:	P = 12.13;	Cc = 0.792892311014899;	K0c = 0.365728618362282;	break;
			//case 4:	P = 18.06;	Cc = 0.670143488818985;	K0c = 0.408018542395396;	break;
			//case 5:	P = 24.10;	Cc = 0.561589528283278;	K0c = 0.458322835247628;	break;
			//case 6:	P = 29.09;	Cc = 0.482872449849904;	K0c = 0.506776012431355;	break;
			case 1:	P =  7.27;	Cc = 0.907367789194129;	K0c = 0.335274865837589;	break;
			case 2:	P = 12.13;	Cc = 0.792892311014899;	K0c = 0.365728618362282;	break;
			case 3:	P = 18.06;	Cc = 0.670143488818985;	K0c = 0.408018542395396;	break;
			case 4:	P = 24.10;	Cc = 0.561589528283278;	K0c = 0.458322835247628;	break;
			case 5:	P = 29.09;	Cc = 0.482872449849904;	K0c = 0.506776012431355;	break;

			/* using normal settings (in vlt_K0cFind.c, vlt_ThermStates.c, and fits) //
			case 1:	P =  0.050;	Cc = 1.10500797392886;	K0c = 0.295614012972922;	break;
			case 2:	P =  1.646;	Cc = 1.0581889918361;	K0c = 0.303871530617906;	break;
			case 3:	P =  7.328;	Cc = 0.905918431108949;	K0c = 0.335619263473799;	break;
			case 4:	P = 15.031;	Cc = 0.730685115724445;	K0c = 0.385639221938482;	break;
			case 5:	P = 18.180;	Cc = 0.667831302441093;	K0c = 0.408941994995917;	break;
			case 6:	P = 22.533;	Cc = 0.588281338805611;	K0c = 0.44447352661157;		break;
			case 7:	P = 25.868;	Cc = 0.532635362503492;	K0c = 0.474695910004971;	break;*/
			/* using B' = 2B (in vlt_K0cFind.c, vlt_ThermStates.c, and fits) //
			case 1:	P =  0.050;	Cc = 1.1003376541326;	K0c = 0.338730088317275;	break;
			case 2:	P =  1.646;	Cc = 1.05648724043934;	K0c = 0.348372781936236;	break;
			case 3:	P =  7.328;	Cc = 0.913079877533694;	K0c = 0.385556454984416;	break;
			case 4:	P = 15.031;	Cc = 0.746257959253795;	K0c = 0.444486928182064;	break;
			case 5:	P = 18.180;	Cc = 0.685924754942954;	K0c = 0.472045003545693;	break;
			case 6:	P = 22.533;	Cc = 0.609090754516276;	K0c = 0.51419404658953;		break;
			case 7:	P = 25.868;	Cc = 0.554925448693439;	K0c = 0.550204504048639;	break; */
		}
		Tc = Tc0 + Tc1*P + Tc2*P*P + Tc3*P*P*P + Tc4*P*P*P*P;

		// Prepare output file, print identification, values, and headings for data (to screen, too) //
		if(TempLabelOp==0)  asprintf(&filename, "3Dvlt_macro_cp_P%06.3f_DK0%1.0e_T_Op%i.dat",  P, DK0, Op);
		if(TempLabelOp==1)  asprintf(&filename, "3Dvlt_macro_cp_P%06.3f_DK0%1.0e_Tv_Op%i.dat", P, DK0, Op);
		outfile = fopen(filename,"w");  //  E.g., "vlt_HeatCap_P.out" or "vlt_HeatCap_P_00.050_DK0_1e-07.dat"
		fprintf(outfile,"# Filename: %s\n", filename);
		fprintf(outfile,"# Source: 3Dvlt_macro_cp.c\n");
		fprintf(outfile,"# Source version: %s\n", "9 (2011 Dec 14 - ...)");
		fprintf(outfile,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g, B=%21.21g, THETA=%g, R=%g\n", PI,PISQ,PICU,B,THETA,R);
		fprintf(outfile,"# Parameter values: N=%i, Kthresh=%g, DK0=%g, dl=%g\n", N,Kthresh,DK0,dl);
		fprintf(outfile,"# Parameter values: P=%g, Cc=%1.15g, K0c=%1.15g, Tc=%g\n", P,Cc,K0c,Tc);
		fprintf(outfile,"# Option values: Op=%i\n", Op);
		fprintf(outfile,"# %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","Kr","K","y","e","amax","Kr/K0","K/K0","G","cp");
		printf(         "\n%s%g\t%s%g\t%s%g\t%s%g%s\n\n", " ---  P = ",P,"Cc = ",Cc,"K0c = ",K0c,"Tc = ",Tc,"  ---");
		printf(         "# %s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","e","Kr/K0","cp");

		// Perform calculations for various temperatures (calculations for each temp are independent) //
		i=0;
		notdone = true;
		while(notdone==true){
			switch(Op){
				case 1:	// A few hand-picked temperatures //
					switch(i){
						case 0:   Tfrac=1.00;  Tempv=1.0-Tfrac;  i++;            break;
						case 1:   Tfrac=0.90;  Tempv=1.0-Tfrac;  i++;            break;
						case 2:   Tfrac=0.80;  Tempv=1.0-Tfrac;  i++;            break;
						case 3:   Tfrac=0.70;  Tempv=1.0-Tfrac;  i++;            break;
						case 4:   Tfrac=0.60;  Tempv=1.0-Tfrac;  i++;            break;
						case 5:   Tfrac=0.50;  Tempv=1.0-Tfrac;  i++;            break;
						case 6:   Tfrac=0.40;  Tempv=1.0-Tfrac;  i++;            break;
						case 7:   Tfrac=0.30;  Tempv=1.0-Tfrac;  i++;            break;
						case 8:   Tfrac=0.20;  Tempv=1.0-Tfrac;  i++;            break;
						case 9:   Tfrac=0.10;  Tempv=1.0-Tfrac;  notdone=false;  break;
					}
					break;
				case 2:	// Decrement below Tc (staying very close to Tc) //
					// from T/Tc = 1-1e-8 to T/Tc = 0.992
					Tempv = 0.0+pow(10,-8.0+0.1*i);  Tfrac = 1.0-Tempv;
					if(i==59) notdone=false;
					i++;  break;
				case 3:	// Increment from near zero to above Tc //
					// from T/Tc = 0.206 to T/Tc = 0.9990;  then from T/Tc = 1.0012 to T/Tc = 2
					if(i<30){ Tempv = 0.0+pow(10,-0.1-0.1*i);      Tfrac = 1.0-Tempv; }
					else    { Tempv = 0.0-pow(10,-5.0+0.1*(i-9));  Tfrac = 1.0-Tempv; }
					if(i==59) notdone=false;
					i++;  break;
				case 4:	// Increment from near zero to just below Tc //
					// from T/Tc = 0.206 to T/Tc = 0.9991
					Tempv = 0.0+pow(10,-0.1-0.05*i);  Tfrac = 1.0-Tempv;
					if(i==59) notdone=false;
					i++;  break;
				case 5:	// Decrement below Tc (staying extremely close to Tc) //
					// T/Tc = 1-1e-20  to  T/Tc ~ 1-1e-16
					Tempv = 0.0+pow(10,-20.0+0.1*i);  Tfrac = 1.0-Tempv;
					// earlier: T/Tc = 1-1e-12  to  T/Tc ~ 1-1e-6
					if(i==59) notdone=false;
					i++;  break;
				case 6:	// Decrement below Tc (staying very close to Tc) //
					// from T/Tc = 1-1e-10 to T/Tc = 0.99992
					Tempv = 0.0+pow(10,-10.0+0.1*i);  Tfrac = 1.0-Tempv;
					if(i==59) notdone=false;
					i++;  break;
				case 7:	// Increment from 0 to 1.1*Tc //
					// from T/Tc = 0 to T/Tc = 2
					Tfrac = 0.001*i;  Tempv = 1.0-Tfrac;
					if(i==2000) notdone=false;
					i++;  break;
			}

			// Set temperature- and pressure-dependent values //
			estimated.tempv = Tempv;
			T = Tc*(1-estimated.tempv);
			rho = rho00         + rho10*T         + rho20*T*T         + rho30*T*T*T         + rho40*T*T*T*T 
			    + rho01*P       + rho11*T*P       + rho21*T*T*P       + rho31*T*T*T*P       + rho41*T*T*T*T*P 
			    + rho02*P*P     + rho12*T*P*P     + rho22*T*T*P*P     + rho32*T*T*T*P*P     + rho42*T*T*T*T*P*P 
			    + rho03*P*P*P   + rho13*T*P*P*P   + rho23*T*T*P*P*P   + rho33*T*T*T*P*P*P   + rho43*T*T*T*T*P*P*P 
			    + rho04*P*P*P*P + rho14*T*P*P*P*P + rho24*T*T*P*P*P*P + rho34*T*T*T*P*P*P*P + rho44*T*T*T*T*P*P*P*P;
			DrhoDT = rho10         + 2.0*rho20*T         + 3.0*rho30*T*T         + 4.0*rho40*T*T*T 
			       + rho11*P       + 2.0*rho21*T*P       + 3.0*rho31*T*T*P       + 4.0*rho41*T*T*T*P 
			       + rho12*P*P     + 2.0*rho22*T*P*P     + 3.0*rho32*T*T*P*P     + 4.0*rho42*T*T*T*P*P 
			       + rho13*P*P*P   + 2.0*rho23*T*P*P*P   + 3.0*rho33*T*T*P*P*P   + 4.0*rho43*T*T*T*P*P*P 
			       + rho14*P*P*P*P + 2.0*rho24*T*P*P*P*P + 3.0*rho34*T*T*P*P*P*P + 4.0*rho44*T*T*T*P*P*P*P;
			D2rhoDT2 = 2.0*rho20         + 6.0*rho30*T         + 12.0*rho40*T*T 
			         + 2.0*rho21*P       + 6.0*rho31*T*P       + 12.0*rho41*T*T*P 
			         + 2.0*rho22*P*P     + 6.0*rho32*T*P*P     + 12.0*rho42*T*T*P*P 
			         + 2.0*rho23*P*P*P   + 6.0*rho33*T*P*P*P   + 12.0*rho43*T*T*P*P*P 
			         + 2.0*rho24*P*P*P*P + 6.0*rho34*T*P*P*P*P + 12.0*rho44*T*T*P*P*P*P;
			a0 = 549.0*Tc*K0c/rho;  // in angstroms
			// a0(T,P) = ( m^2 kB Tc(P) K0c(P) ) / (hbar^2 rho(T,P) )
			// where m = 6.65e-27 is the mass of 4He atom in kg, kB is the Boltzmann constant, and hbar is the Dirac constant
			// COULD UPDATE THAT:  use  m = (standard atomic weight)/(Avogadro constant) = (0.004002602 kg/mol)/(6.02214179 x 10^23 atoms/mol) = 6.64647585456 x 10^-27
			//                     so   m^2 kB / hbar^2 = (6.64647585456 x 10^-27)^2 * (1.3806503 x 10^-23) / ( (6.626068 x 10^-34)/2pi )^2 = 5.48423431589 x 10^-8
			//                     and  a0 = (548.423431589)*Tc*K0c/rho;  // in angstroms
			// AND UPDATE eqns below, including cP
			Da0DT = -549.0*Tc*K0c*DrhoDT/(rho*rho);
			D2a0DT2 = 2.0*549.0*Tc*K0c*DrhoDT*DrhoDT/(rho*rho*rho) - 549.0*Tc*K0c*D2rhoDT2/(rho*rho);

			// Initialize lengthscale, energy values, and in/decrements at zero //
			estimated.l = 0.0;
			plus.l      = 0.0;
			minus.l     = 0.0;
			estimated.thrmv[2] = 0.0;
			plus.thrmv[2]      = 0.0;
			minus.thrmv[2]     = 0.0;

			// Initialize temperature-dependent values at initial lengthscale //
			K0 = K0c/(1.0-estimated.tempv);
			estimated.thrmv[0] = K0;
			plus.thrmv[0]      = K0+DK0;
			minus.thrmv[0]     = K0-DK0;
			estimated.thrmv[1] = exp(-PISQ*K0*Cc);
			plus.thrmv[1]      = exp(-(K0+DK0)*PISQ*Cc);
			minus.thrmv[1]     = exp(-(K0-DK0)*PISQ*Cc);

			// Calculate state-variables at progressively larger length scales, up to coherence length, via Runge-Kutta method (derk) //
			//while(estimated.thrmv[0]<Kthresh){
			while(estimated.thrmv[0]<Kthresh && estimated.l<=lmax){
				// Calculate next estimated state and different temperature (K0-value) test states //
				// (to enable a derivative with respect to K0), incrementing lengthscale           //
				estimated = derk(IN vltRecRel, INOUT estimated, IN dl, IN N);
				plus      = derk(IN vltRecRel, INOUT plus,      IN dl, IN N);
				minus     = derk(IN vltRecRel, INOUT minus,     IN dl, IN N);
				/* Print intermediate numbers to screen (for bug checking)
				j++;
				if(j%2000==0){
					Kr = estimated.thrmv[0]*exp(-estimated.l);
					DeDK0plus = (plus.thrmv[2]-estimated.thrmv[2])/DK0;
					DeDK0minus = (estimated.thrmv[2]-minus.thrmv[2])/DK0;
					DeDK0 = (DeDK0plus+DeDK0minus)/2;
					D2eDK02 = (DeDK0plus-DeDK0minus)/DK0;
					cp = -R*estimated.thrmv[2]*((6.65e3)/(rho*a0*a0*a0))*T*((T*alpha-2)*alpha-T*DalphaDT+2*(T*alpha-1)*Da0DT/a0+T*(Da0DT*Da0DT)/(a0*a0)-T*D2a0DT2/a0)
					  - R*((6.65e3)/(rho*a0*a0*a0))*((T*alpha+2*T*Da0DT/a0)*K0*DeDK0+K0*K0*D2eDK02)
					  - T*P*(4.005e-4/rho)*(alpha*alpha-DalphaDT);
					printf(         "%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\n",
						1-estimated.tempv, estimated.tempv, estimated.l, estimated.thrmv[0], estimated.thrmv[1], estimated.thrmv[2], Kr, Kr/K0, cp);
				}
				if(j>10000) j=0; */
			}

			// Update/set values that depend on state-variables, if not using "Print intermediate numbers..." block of code //
			amax = a0*exp(estimated.l);
			Kr   = estimated.thrmv[0]*exp(-estimated.l);
			G    = (4.0/PISQ)*estimated.thrmv[1]/(amax*amax*amax*amax*amax*amax);
			DeDK0plus  = (plus.thrmv[2]-estimated.thrmv[2])/DK0;
			DeDK0minus = (estimated.thrmv[2]-minus.thrmv[2])/DK0;
			DeDK0   = (DeDK0plus+DeDK0minus)/2;
			D2eDK02 = (DeDK0plus-DeDK0minus)/DK0;

			// Calculate individual terms of cp (in terms of rho) //
			cp01 = R*((6.65e3)/(rho*a0*a0*a0))*K0*K0*D2eDK02*( -1.0 );
			cp02 = R*((6.65e3)/(rho*a0*a0*a0))*K0*DeDK0*(  4.0*(T/rho)*DrhoDT      );
			cp03 = R*((6.65e3)/(rho*a0*a0*a0))*estimated.thrmv[2]*( -8.0*(T/rho)*DrhoDT                );
			cp04 = R*((6.65e3)/(rho*a0*a0*a0))*estimated.thrmv[2]*( -2.0*(T*T/(rho*rho))*DrhoDT*DrhoDT );
			cp05 = R*((6.65e3)/(rho*a0*a0*a0))*estimated.thrmv[2]*( -2.0*(T*T/rho)*D2rhoDT2            );
			cp06 = -T*(100000*P)*(4.005e-4/rho)*(  2.0*(1.0/(rho*rho))*DrhoDT*DrhoDT );
			cp07 = -T*(100000*P)*(4.005e-4/rho)*( -(1.0/rho)*D2rhoDT2                );

			// Calculate cp //
			cp = cp01+cp02+cp03+cp04+cp05+cp06+cp07;

			// Print data to screen and output file //
			fprintf(outfile,"%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\n",
				1-estimated.tempv, estimated.tempv, estimated.l, Kr, estimated.thrmv[0], estimated.thrmv[1], estimated.thrmv[2], amax, Kr/K0, estimated.thrmv[0]/K0, G, cp);
			printf(         "%5.4g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\t%17.15g\n",
				1-estimated.tempv, estimated.tempv, estimated.l, estimated.thrmv[2], Kr/K0, cp);
		}

		// Close output file //
		fclose(outfile);
	}
}



// vortex loop theory recursion relations //
RETARRAY vltRecRel(double l, double z[]){
	RETARRAY dzdl;

	dzdl.arr[0] = z[0]-B*z[0]*z[0]*z[1];
	dzdl.arr[1] = z[1]*( 6.0 - PISQ*z[0]*(1.0-THETA*log(z[0])) );
	dzdl.arr[2] = -PI*z[1]*exp(-3.0*l);

	return dzdl;
}



//  Program Notes  //
//=========================================================================//
/*

See further documentation (UCLAresearchNotes.pdf, UCLAresearchPrograms.pdf, derk1.c) for elaboration.


plus and minus are test states (with test state-variables)

*/
