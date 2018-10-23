//  File Comments  //
//=========================================================================//

/* FILENAME: 3Dvlt_constT.c */
/* VERSION: 10 (2012 Jul 29 - ...)
     This program was formerly called vlt_ThermStates.c.
     Got rid of Kthresh cutoff.
     Added l=0 data print out.  */
/* AUTHOR: Andrew Forrester <aforrester@ucla.edu> */
/* DESCRIPTION:
   * This C program calculates (length-scale-dependent) thermodynamic states of a system of liquid helium-4 (4He) under conditions specified by user inputs.
   * The system is a spherical region (of variable diameter greater than D) of liquid 4He at some constant pressure, near or below the critical transition temperature Tc (which depends on pressure) separating normal fluid and superfluid (aka the lambda temperature).
   * A state is merely a set of values of certain quantities that correspond to temperature T, length scale a, dimensionless superfluid ratio K, dimensionless fugacity y, and dimensionless Helmholtz energy parameter e.  The "bare" superfluid ratio K0 at the smallest relevant length scale a0 (the diameter of the smallest vortex loops) and the renormalized superfluid ratio Ks (and Ks/K0) are also of interest.  See the "Program Notes" below the code for explanation of these quantities and the constants, variables, and functions in this program.
   * The actual values of the set of temperatures {Tn} and the length scale a are irrelevant; the calculations are made using unitless variables (tempv and l) scaled by Tc (which depends on pressure) and a0 (which depends on temperature and pressure).  The actual value of the pressure is also irrelevant; the user may choose a value of Cc without knowing what pressure it corresponds to.
   * The program uses Gary A. Williams's vortex-loop theory of the superfluid phase transition (i.e., its differential recursion relations) and a 4th-order Runge-Kutta integration method to calculate (length-scale-dependent) thermodynamic quantities (essentially K0,K,Ks,y,e) of the system at a set of temperatures {Tn}, at various length scales (l, no greater than a length scale corresponding to the choice of lmax), and at some pressure (corresponding to the choice of Cc, the critical dimensionless vortex-loop-core energy).
   * For a set of temperatures {Tn}, the program starts from some initial state (l=0,K0=K0(Cc,Tn),K=K0,Ks=K0,y=y(K0,Cc),e=0) (depending on each temperature Tn and the choice of Cc) at a smallest length scale (l=0 and a=a0), it integrates to larger length scales (incrementing by Dl using the Runge-Kutta technique) up to some limiting length scale (see while loop, where l<lmax and K<5), and it estimates and records the state at that final length scale (l=lf,K0=K0(Cc,Tn),K=Kf,Ks=Ks(K,l),y=yf,e=ef).
   * The various temperatures are determined by "Options" (AA, A, B, C, and D) in the program, where one can choose (by manually editing this file) sets of temperatures below, above, or very near Tc.  One could also add more options to choose different sets of temperatures.

   * In the basename "vlt_ThermStates", "vlt" refers to Williams's Vortex Loop Theory (of the superfluid phase transition) and "ThermStates" refers to the THERModynamic STATES described by the calculated thermodynamic quantities.  (This file was previously called ring3.c.)

   Inputs:
   * N       - You can choose N=2 (if you don't care to calculate e) or N=3.
   * Kthresh -
   * dl      -
   * lmax    - Choosing 100 corresponds to a system with diameter no greater than D = a0*exp(lmax) ~ 6.7e33 meters (given that a0 ~ 2.5 angstroms at low pressures), allowing for a huge system and calculations at enormous length scales.
   * From vlt_K0cFind.c, you can get a set of {Cc, K0c} data to put in the 
   * Select (uncomment) "Option" (AA, A, B, C, or D) to choose a set of temperatures.

   Output:
   * This program returns text files of data in a two dimensional arrays.

   Questions:
   1. Should we say that the system is of variable volume or should we say that the system is of infinite volume?  (D would then just be a limiting diameter for the calculations.)  If we are assuming that the system is of fixed, finite volume (diameter D), then when you select the temperature of the system, doesn't that also determine the pressure (and thus, wouldn't that determine a certain Cc for each temperature)?
   2. Why do you use
      estimated.thrmv[0]<5.0
   in this program but
      estimated.thrmv[0]<6.0
   in vlt_HeatCap_P.c?
   3. Why don't you include a finite-size option in vlt_HeatCap_P.c
      while(estimated.thrmv[0]<6)
   like you do in this file:
      while(estimated.l<lmax && estimated.thrmv[0]<5.0) ?
*/
/* EXT FILES: none */
/* COMPILE NOTES:
   * To compile, type "g++ -lm 3Dvlt_constT.c" without the quotes; then to run, type "./a.out".
   * WARNING: Make sure that the number of files generated is not unexpectedly large. (#files = #CcK0cSets X #temps)
   * NOTE: Be sure lsteps is divisible by ldatamax.
*/
/* TROUBLESHOOTING IDEAS:
   * Take smaller step sizes (lsteps and time steps dt0)
   * Maybe change Kthresh to something greater than 5.0 ?
*/



//  Function Preparation  //
//=========================================================================//

// Standard routine header files //
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Constants //
const double PI    = 3.14159265358979323846;
const double PISQ  = 9.86960440108935861883;
const double PICU  = 31.00627668029981620634;
const double A3    = 4.0*PICU/3.0;
const double THETA = 0.6;

// Program parameters/inputs //
const int    OpA  = 1;   //  Option A   Set recursion relation behavior
                         //              1 -> "itlog flow": use if-then statement to limit ac at a -- "if-then log flow" of EqRecRel1;
                         //              2 -> "w/log flow": no if-then, using logarithm so ac may exceed a (two-fixed points) -- "with log flow" of EqRecRel2;
                         //              3 -> "nolog flow": no if-then, no logarithm so ac = a the whole time (one fixed pt) -- "no log flow" of EqRecRel3;
                         //              4 -> "w/log flow" using Kthresh: no if-then, using logarithm so ac may exceed a (two-fixed points) -- "with log flow" of EqRecRel2.
const int    OpB  = 1;   //  Option B   Temperatures option
                         //              1 -> hand-picked temp's; else temperature-spread formulas.
const int    N        = 3;         //  Selecting quantities to be calculated (2 -> K and y/G; 3 -> K, y/G, and e)
const double lmax     = 100.0;      //  1  //  maximum of the length scale l (potentially approximate)
const int    lsteps   = 10000000;   //  1250, 5000, 100000 //  lsteps = lmax/dl = 10/0.0001 = 100,000  (make sure this is divisible by ldatamax)
const int    lpts     = lsteps+1;  //  from l=0 to l=lmax, inclusive (incl. the boundary-condition-enforcing pnt, there are lpts+1 "l" pnts)
const int    ldatamax = 1000;      //  50,100 //  max number of data points recorded per temperature examined
const double a0       = 1.0;       //  a0 in units of a0 (!)
const double a03      = a0*a0*a0;  //  a0 to the third power
const double a06      = a03*a03;   //  a0 to the sixth power
const double KstarA   = 4.146766416946759293;
const double Kthresh  = KstarA;       //  (used to use value of 5.0 or 6.0; use this when NOT using the if-then code in EqRecRel)



//  Function Prototypes  //
//=========================================================================//

// rk4((*f)(),*x,y[],h,n);
// EqRecRel(x,z,dzdx,n);
void rk4(void (*f)(double, double*, double*, unsigned int), double *x, double y[], double h, unsigned n);
void EqRecRel1(double x, double z[], double dzdx[], unsigned n);
void EqRecRel2(double x, double z[], double dzdx[], unsigned n);
void EqRecRel3(double x, double z[], double dzdx[], unsigned n);



//  Function Definitions  //
//=========================================================================//

int main(void){
	// Main function definitions //
	int    s,i,j, notdone=true, TempLabelOp;
	int    lstep, datapt;
	double dblsteps=lsteps, dl=lmax/dblsteps;
	int    ldataskip = lsteps/ldatamax;
	double l, r, rc; // r=a/a0, rc = ac/a
	double Tfrac, tempv;
	double z[N], K, y, e, Ks, G;
	double K0, y0;
	double Cc, K0c;
	FILE   *outfile;
	char   *filename;
	void   (*EqRecRel)(double x, double z[], double dzdx[], unsigned n);

	switch(OpA){ // Set recursion relation behavior
		case 1:  // "itlog flow": use if-then statement to limit ac at a -- "if-then log flow"
			 EqRecRel = EqRecRel1;
			 break;
		case 2:  // "w/log flow": no if-then, using logarithm so ac may exceed a (two-fixed points) -- "with log flow"
			 EqRecRel = EqRecRel2;
			 break;
		case 3:  // "nolog flow": no if-then, no logarithm so ac = a the whole time (one fixed pt) -- "no log flow"
			 EqRecRel = EqRecRel3;
			 break;
		case 4:  // "w/log flow" using Kthresh: no if-then, using logarithm so ac may exceed a (two-fixed points) -- "with log flow" of EqRecRel2;
			 EqRecRel = EqRecRel2;
			 break;
	}

	// Progress through sets of values of Cc and K0c, where a data file will be produced for each set //
	//for(s=1; s<=21; s++){
	//for(s=1; s<=22; s++){
	//for(s=1; s<=10; s++){
	//for(s=1; s<=1; s++){
	for(s=1; s<=2; s++){
		switch(s){
			// from plots_earlyequil/plots09/fit_Pdep_Cc_results.dat and 3Dvlt_K0cFind.c //
			//case 1:	Cc = 1.10500797392886;	K0c = 0.295614012972922;	break;  // P = 0.050;

			/* Just checking for ac/a = rc */
			case 1:		Cc = 0.40;	K0c = 0.57534342964893814;	break;  // ac/a = K^theta*exp(C) = 1.070715 (@ Tc, l=0)
			case 2:		Cc = 1.10;	K0c = 0.29646893900678273;	break;  // ac/a = K^theta*exp(C) = 1.448476 (@ Tc, l=0)

			/* Best yet calculation (using "itlog flow"), checking to see if it matches thesis values (which used Kthresh)
			// And (it turns out) I'll have to see how the Kthresh values (using this new code) turn out in comparison, since the results came out funny
			// from 3Dvlt_K0cFind.c v11 - 3Dvlt_K0cFind_Output38_uncert1e-17_dl1e-05.dat - "itlog flow" (if-then statement used in recursion relations)
			case 1:		Cc = 0.20;	K0c = 0.9346549781218052;	break;  // ac/a = K^theta*exp(C) = 1.172869 (@ Tc, l=0)
			case 2:		Cc = 0.25;	K0c = 0.79674794984036801;	break;  // ac/a = K^theta*exp(C) = 1.120382 (@ Tc, l=0)
			case 3:		Cc = 0.30;	K0c = 0.70109094206651523;	break;  // ac/a = K^theta*exp(C) = 1.090820 (@ Tc, l=0)
			case 4:		Cc = 0.35;	K0c = 0.63025269735391265;	break;  // ac/a = K^theta*exp(C) = 1.075751 (@ Tc, l=0)
			case 5:		Cc = 0.40;	K0c = 0.57534342964893814;	break;  // ac/a = K^theta*exp(C) = 1.070715 (@ Tc, l=0)
			case 6:		Cc = 0.45;	K0c = 0.53132345421063487;	break;  // ac/a = K^theta*exp(C) = 1.073119 (@ Tc, l=0)
			case 7:		Cc = 0.50;	K0c = 0.49510828997129153;	break;  // ac/a = K^theta*exp(C) = 1.081353 (@ Tc, l=0)
			case 8:		Cc = 0.55;	K0c = 0.46469704703948356;	break;  // ac/a = K^theta*exp(C) = 1.094369 (@ Tc, l=0)
			case 9:		Cc = 0.60;	K0c = 0.43873100720691904;	break;  // ac/a = K^theta*exp(C) = 1.111465 (@ Tc, l=0)
			case 10:	Cc = 0.65;	K0c = 0.41625279991888953;	break;  // ac/a = K^theta*exp(C) = 1.132155 (@ Tc, l=0)
			case 11:	Cc = 0.70;	K0c = 0.39656697289619258;	break;  // ac/a = K^theta*exp(C) = 1.156102 (@ Tc, l=0)
			case 12:	Cc = 0.75;	K0c = 0.37915524531909928;	break;  // ac/a = K^theta*exp(C) = 1.183072 (@ Tc, l=0)
			case 13:	Cc = 0.80;	K0c = 0.36362286558404033;	break;  // ac/a = K^theta*exp(C) = 1.212904 (@ Tc, l=0)
			case 14:	Cc = 0.85;	K0c = 0.3496634698826776;	break;  // ac/a = K^theta*exp(C) = 1.245491 (@ Tc, l=0)
			case 15:	Cc = 0.90;	K0c = 0.33703537107020987;	break;  // ac/a = K^theta*exp(C) = 1.280768 (@ Tc, l=0)
			case 16:	Cc = 0.95;	K0c = 0.32554514449601912;	break;  // ac/a = K^theta*exp(C) = 1.318701 (@ Tc, l=0)
			case 17:	Cc = 1.00;	K0c = 0.31503600689857347;	break;  // ac/a = K^theta*exp(C) = 1.359285 (@ Tc, l=0)
			case 18:	Cc = 1.05;	K0c = 0.3053794236081705;	break;  // ac/a = K^theta*exp(C) = 1.402533 (@ Tc, l=0)
			case 19:	Cc = 1.10;	K0c = 0.29646893900678273;	break;  // ac/a = K^theta*exp(C) = 1.448476 (@ Tc, l=0)
			case 20:	Cc = 1.15;	K0c = 0.28821556877326282;	break;  // ac/a = K^theta*exp(C) = 1.497163 (@ Tc, l=0)
			case 21:	Cc = 1.20;	K0c = 0.28054430897522864;	break;  // ac/a = K^theta*exp(C) = 1.548653 (@ Tc, l=0) */

			/* Testing the renormalization flow idea: if I set lmax very high, will I see deviations in the expected pattern?
			// Also, does the 3Dvlt_K0cFind3_Plot.c method of finding K0c change the result much? //
			// Values from 3Dvlt_K0cFind.c
			case 1:	Cc = 1.10;	K0c = 0.296468939006784;	break;
			case 2:	Cc = 0.40;	K0c = 0.575343429658511;	break;
			// Values from 3Dvlt_K0cFind3_Plot.c
			//3Dvlt_K0cFind3_Plot_Output2_OpC2_uncert1e-10.dat  const dl
			//case 1:	Cc = 1.10;	K0c = 0.29646894595846918;	break;
			//case 2:	Cc = 0.40;	K0c = 0.57534343604480565;	break;
			//3Dvlt_K0cFind3_Plot_Output7_OpC4_uncert1e-15.dat  limtd ds
			//case 1:	Cc = 1.10;	K0c = 0.29646893900647414;	break;
			//case 2:	Cc = 0.40;	K0c = 0.5753434296488823;	break; */

			/* from vlt_K0cFind.c, normal settings //
			case 1:		Cc = 1.20;	K0c = 0.280544308975229;	break;
			case 2:		Cc = 1.10;	K0c = 0.296468939006784;	break;
			case 3:		Cc = 1.06;	K0c = 0.30354071294185;		break;
			case 4:		Cc = 1.05;	K0c = 0.305379423608169;	break;
			case 5:		Cc = 1.04;	K0c = 0.307247922776542;	break;
			case 6:		Cc = 1.03;	K0c = 0.309146998428258;	break;
			case 7:		Cc = 1.02;	K0c = 0.311077467512369;	break;
			case 8:		Cc = 1.01;	K0c = 0.313040177313466;	break;
			case 9:		Cc = 1.00;	K0c = 0.315036006898571;	break;
			case 10:	Cc = 0.99;	K0c = 0.317065868648624;	break;
			case 11:	Cc = 0.98;	K0c = 0.319130709880376;	break;
			case 12:	Cc = 0.97;	K0c = 0.321231514565087;	break;
			case 13:	Cc = 0.90;	K0c = 0.337035371070207;	break;
			case 14:	Cc = 0.80;	K0c = 0.36362286558404;		break;
			case 15:	Cc = 0.70;	K0c = 0.39656697289679;		break;
			case 16:	Cc = 0.60;	K0c = 0.438731007210038;	break;
			case 17:	Cc = 0.55;	K0c = 0.464697047043989;	break;
			case 18:	Cc = 0.50;	K0c = 0.495108289977297;	break;
			case 19:	Cc = 0.40;	K0c = 0.575343429658511;	break;
			case 20:	Cc = 0.30;	K0c = 0.701090942081052;	break;
			case 21:	Cc = 0.20;	K0c = 0.934654978144903;	break;
			case 22:	Cc = 0.10;	K0c = 1.52472890244101;		break;*/
			/* from vlt_K0cFind.c using vltRecRel with A3' = 2A3: // redo them */
			/* from vlt_K0cFind.c using vltRecRel with A3' = A3/2: // redo them */
		}
		printf(         "\n%s%g\t%s%g%s\n\n", " ---  Cc = ",Cc,"K0c = ",K0c,"  ---");

		// Select temperatures (calculations for each temp are independent), where a data file will be produced for each temp //
		i=0;
		notdone=true;
		while(notdone==true){
			switch(OpB){
				case 1:	// A few hand-picked temperatures //
					switch(i){
						/*
						case 0:   Tfrac=1.010;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 1:   Tfrac=1.005;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 2:   Tfrac=1.000;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 3:   Tfrac=0.500;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 4:   Tfrac=0.100;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 5:   Tfrac=0.050;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 6:   Tfrac=0.020;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 7:   Tfrac=0.017;  tempv=1.0-Tfrac;  notdone=false;  TempLabelOp=0;  break; */
						/*
						case 0:  Tfrac=1.00;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 1:  Tfrac=0.90;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 2:  Tfrac=0.80;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 3:  Tfrac=0.70;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 4:  Tfrac=0.60;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 5:  Tfrac=0.50;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 6:  Tfrac=0.40;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 7:  Tfrac=0.30;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 8:  Tfrac=0.20;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 9:  Tfrac=0.10;  tempv=1.0-Tfrac;  notdone=false;  TempLabelOp=0;  break; */
						/*
						case 0:  Tfrac=1.00;    tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 1:  Tfrac=0.9999;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 2:  Tfrac=0.999;   tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 3:  Tfrac=0.99;    tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 4:  Tfrac=0.90;    tempv=1.0-Tfrac;  notdone=false;  TempLabelOp=0;  break; */
						/* */
						case 0:  tempv=0.0;       Tfrac=1.0-tempv;  i++;            TempLabelOp=1;  break;
						case 1:  tempv=25.0e-6;   Tfrac=1.0-tempv;  i++;            TempLabelOp=1;  break;
						case 2:  tempv=50.0e-6;   Tfrac=1.0-tempv;  i++;            TempLabelOp=1;  break;
						case 3:  tempv=75.0e-6;   Tfrac=1.0-tempv;  i++;            TempLabelOp=1;  break;
						case 4:  tempv=100.0e-6;  Tfrac=1.0-tempv;  notdone=false;  TempLabelOp=1;  break;
					}
					break;
				// Before using the below code... DO I REALLY WANT 60 FILES??!!
				case 2:	// Decrement below Tc (staying very close to Tc) //
					// from T/Tc = 1-1e-8 to T/Tc = 0.992
					tempv = 0.0+pow(10,-8.0+0.1*i);
					if(i==59) notdone=false;
					i++;  TempLabelOp=1;  break;
				case 3:	// Increment from near zero to above Tc //
					// from T/Tc = 0.206 to T/Tc = 0.9990;  then from T/Tc = 1.0012 to T/Tc = 2
					if(i<30) tempv = 0.0+pow(10,-0.1-0.1*i);
					else     tempv = 0.0-pow(10,-5.0+0.1*(i-9));
					if(i==59) notdone=false;
					i++;  TempLabelOp=0;  break;
				case 4:	// Increment from near zero to just below Tc //
					// from T/Tc = 0.206 to T/Tc = 0.9991
					tempv = 0.0+pow(10,-0.1-0.05*i);
					if(i==59) notdone=false;
					i++;  TempLabelOp=0;  break;
				case 5:	// Decrement below Tc (staying extremely close to Tc) //
					// T/Tc = 1-1e-20  to  T/Tc ~ 1-1e-16
					tempv = 0.0+pow(10,-20.0+0.1*i);
					// earlier: T/Tc = 1-1e-12  to  T/Tc ~ 1-1e-6
					if(i==59) notdone=false;
					i++;  TempLabelOp=1;  break;
				case 6:
					switch(i){
						case 0:  tempv=pow(10,-0.1);  Tfrac=1.0-tempv;  i++;           TempLabelOp=0;  break; // Tfrac=0.206
						case 1:  tempv=pow(10,-1);    Tfrac=1.0-tempv;  i++;           TempLabelOp=0;  break; // Tfrac=0.9
						case 2:  tempv=pow(10,-2);    Tfrac=1.0-tempv;  i++;           TempLabelOp=0;  break; // Tfrac=0.99
						case 3:  tempv=pow(10,-3);    Tfrac=1.0-tempv;  i++;           TempLabelOp=0;  break; // Tfrac=0.999
						case 4:  tempv=-pow(10,-2.9); Tfrac=1.0-tempv;  i++;           TempLabelOp=0;  break; // Tfrac=1.001
						case 5:  tempv=-pow(10,-2);   Tfrac=1.0-tempv;  i++;           TempLabelOp=0;  break; // Tfrac=1.01
						case 6:  tempv=-pow(10,-1);   Tfrac=1.0-tempv;  i++;           TempLabelOp=0;  break; // Tfrac=1.1
						case 7:  tempv=-pow(10,-0);   Tfrac=1.0-tempv;  notdone=false; TempLabelOp=0;  break; // Tfrac=2
					}
			}

			// Prepare output file(s), print identification and values //
			if(TempLabelOp==0)  asprintf(&filename, "3Dvlt_constT_%4.3f_Op%i%i_Cc%3.2f_lmax%g_dl%g.dat", Tfrac,OpA,OpB,Cc,lmax,dl);
			if(TempLabelOp==1)  asprintf(&filename, "3Dvlt_constTv_%g_Op%i%i_Cc%3.2f_lmax%g_dl%g.dat",   tempv,OpA,OpB,Cc,lmax,dl);
			outfile = fopen(filename,"w");  //  E.g., "3Dvlt_constT.out"
			//printf(     "\n\n# Filename: %s\n", filename);
			printf( "Cc is now %g  and filename = %s\n", Cc,filename);
			fprintf(outfile,"# Filename: %s\n", filename);
			fprintf(outfile,"# Source: 3Dvlt_constT.c\n");
			fprintf(outfile,"# Source version: %s\n", "10 (2012 Jul 29 - ...)");
			fprintf(outfile,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g, A3=%21.21g, THETA=%g\n", PI,PISQ,PICU,A3,THETA);
			fprintf(outfile,"# Parameter values: OpA=%i, OpB=%i, N=%i, Tfrac=%g, tempv=%g, Cc=%g, K0c=%g, lmax=%g, lsteps=%i, (dl=%g), lpts=%i, ldatamax=%i, a0=%g, a03=%g, a06=%g\n", OpA,OpB,N,Tfrac,tempv,Cc,K0c,lmax,lsteps,dl,lpts,ldatamax,a0,a03,a06);
			fprintf(outfile,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "Data-pt","l-Step","T/Tc","1-T/Tc","l","K","y","e","K/K0","Ks/K0","G","a/a0","ac/a");
			printf(         "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "Data-pt","l-Step","T/Tc","1-T/Tc","l","K","y","e","K/K0","Ks/K0","G","a/a0","ac/a");

			// Initialize quantities at smallest length-scale (l=0) //
			l = 0.0;
			lstep = 0;
			z[0]  = K = K0 = K0c/(1.0-tempv);
			z[1]  = y = y0 = exp(-PISQ*K0*Cc);
			z[2]  = e = 0.0;
			Ks = K*exp(-l);
			G  = y*exp(-6.0*l);
			r = exp(l);
			if(K < 1.0) rc = pow(K,THETA);
			else        rc = 1.0; // if K > 1 (ac' = a*K^THETA > a), then set ac' = a
			datapt = 1;
			printf(         "%i\t%i\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
				datapt, lstep, Tfrac, tempv, l, K, y, e, K/K0, Ks/K0, G, r, rc);
			fprintf(outfile,"%i\t%i\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\n",
				datapt, lstep, Tfrac, tempv, l, K, y, e, K/K0, Ks/K0, G, r, rc);

			// Calculate state-variables at progressively larger length scales, up to coherence length, via Runge-Kutta method (rk4) //
			if(OpA!=4){
				while(lstep<lsteps && l<lmax){
					// Step out in length scale l, calculating next l, K, y, and e (and G and Ks) //
					rk4(EqRecRel, &l, z, dl, N);  // equil: rk4 3D (K,y,e), EqRecRel
					lstep++;
					// Bug Check: Print at most [ldatamax] (e.g., 50) data points to screen and output file //
					if(lstep%ldataskip==0){
						K  = z[0];
						y  = z[1];
						e  = z[2];
						Ks = K*exp(-l);
						G  = y*exp(-6.0*l);
						r  = exp(l);
						if(K < 1.0) rc = pow(K,THETA);
						else        rc = 1.0; // if K > 1 (ac' = a*K^THETA > a), then set ac' = a
						datapt++;
						printf(         "%i\t%i\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
							datapt, lstep, Tfrac, tempv, l, K, y, e, K/K0, Ks/K0, G, r, rc);
						fprintf(outfile,"%i\t%i\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\n",
							datapt, lstep, Tfrac, tempv, l, K, y, e, K/K0, Ks/K0, G, r, rc);
					}
				}
			}
			if(OpA==4){
				while(lstep<lsteps && l<lmax && K<Kthresh){
					// Step out in length scale l, calculating next l, K, y, and e (and G and Ks) //
					rk4(EqRecRel, &l, z, dl, N);  // equil: rk4 3D (K,y,e), EqRecRel
					lstep++;
					K  = z[0];
					// Print at most [ldatamax] (e.g., 50) data points to screen and output file //
					if(lstep%ldataskip==0){
						y  = z[1];
						e  = z[2];
						Ks = K*exp(-l);
						G  = y*exp(-6.0*l);
						r  = exp(l);
						if(K < 1.0) rc = pow(K,THETA);
						else        rc = 1.0; // if K > 1 (ac' = a*K^THETA > a), then set ac' = a
						datapt++;
						printf(         "%i\t%i\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
							datapt, lstep, Tfrac, tempv, l, K, y, e, K/K0, Ks/K0, G, r, rc);
						fprintf(outfile,"%i\t%i\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\n",
							datapt, lstep, Tfrac, tempv, l, K, y, e, K/K0, Ks/K0, G, r, rc);
					}
				}
				K  = z[0];
				y  = z[1];
				e  = z[2];
				Ks = K*exp(-l);
				G  = y*exp(-6.0*l);
				r  = exp(l);
				if(K < 1.0) rc = pow(K,THETA);
				else        rc = 1.0; // if K > 1 (ac' = a*K^THETA > a), then set ac' = a

				// Record (final, macroscopic) data values //
				datapt++;
				printf(         "%i\t%i\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
					datapt, lstep, Tfrac, tempv, l, K, y, e, K/K0, Ks/K0, G, r, rc);
				fprintf(outfile,"%i\t%i\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\t%17.17g\n",
					datapt, lstep, Tfrac, tempv, l, K, y, e, K/K0, Ks/K0, G, r, rc);
				fflush(outfile);
			}

			fclose(outfile);
		}
	}

	return EXIT_SUCCESS;
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



// Equilibrium (K,y) recursion relations (version 1, proper) //
void EqRecRel1(double x, double z[], double dzdx[], unsigned n){
	dzdx[0] = z[0]-A3*z[0]*z[0]*z[1];
	if( z[0] > 1.0 ) // if K > 1 (ac' = a*K^THETA > a), then set ac' = a (replace log(K[i]) with zero)
		dzdx[1] = z[1]*( 6.0 - PISQ*z[0] );
	else
		dzdx[1] = z[1]*( 6.0 - PISQ*z[0]*(1.0-THETA*log(z[0])) );
	dzdx[2] = -PI*exp(-3.0*x)*z[1];
}



// Equilibrium (K,y) recursion relations (version 2) //
void EqRecRel2(double x, double z[], double dzdx[], unsigned n){
	dzdx[0] = z[0]-A3*z[0]*z[0]*z[1];
	//if( z[0] > 1.0 ) // if K > 1 (ac' = a*K^THETA > a), then set ac' = a (replace log(K[i]) with zero)
	//	dzdx[1] = z[1]*( 6.0 - PISQ*z[0] );
	//else
		dzdx[1] = z[1]*( 6.0 - PISQ*z[0]*(1.0-THETA*log(z[0])) );
	dzdx[2] = -PI*exp(-3.0*x)*z[1];
}



// Equilibrium (K,y) recursion relations (version 3) //
void EqRecRel3(double x, double z[], double dzdx[], unsigned n){
	dzdx[0] = z[0]-A3*z[0]*z[0]*z[1];
	//if( z[0] > 1.0 ) // if K > 1 (ac' = a*K^THETA > a), then set ac' = a (replace log(K[i]) with zero)
		dzdx[1] = z[1]*( 6.0 - PISQ*z[0] );
	//else
	//	dzdx[1] = z[1]*( 6.0 - PISQ*z[0]*(1.0-THETA*log(z[0])) );
	dzdx[2] = -PI*exp(-3.0*x)*z[1];
}



//  Program Notes  //
//=========================================================================//
/*

...

*/
