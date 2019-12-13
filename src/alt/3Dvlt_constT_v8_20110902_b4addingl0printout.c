//  File Comments  //
//=========================================================================//

/* FILENAME: 3Dvlt_constT.c */
/* VERSION: 8 (2011 Jun 30 - 2011 Sep 02)
     This program was formerly called vlt_ThermStates.c. */
/* AUTHOR: Andrew Forrester <aforrester@ucla.edu> */
/* DESCRIPTION:
   * This C program calculates (length-scale-dependent) thermodynamic states of a system of liquid helium-4 (4He) under conditions specified by user inputs.
   * The system is a spherical region (of variable diameter greater than D) of liquid 4He at some constant pressure, near or below the critical transition temperature Tc (which depends on pressure) separating normal fluid and superfluid (aka the lambda temperature).
   * A state is merely a set of values of certain quantities that correspond to temperature T, length scale a, dimensionless superfluid ratio K, dimensionless fugacity y, and dimensionless Helmholtz energy parameter e.  The "bare" superfluid ratio K0 at the smallest relevant length scale a0 (the diameter of the smallest vortex loops) and the renormalized superfluid ratio Kr (and Kr/K0) are also of interest.  See the "Program Notes" below the code for explanation of these quantities and the constants, variables, and functions in this program.
   * The actual values of the set of temperatures {Tn} and the length scale a are irrelevant; the calculations are made using unitless variables (tempv and l) scaled by Tc (which depends on pressure) and a0 (which depends on temperature and pressure).  The actual value of the pressure is also irrelevant; the user may choose a value of Cc without knowing what pressure it corresponds to.
   * The program uses Gary A. Williams's vortex-loop theory of the superfluid phase transition (i.e., its differential recursion relations) and a 4th-order Runge-Kutta integration method to calculate (length-scale-dependent) thermodynamic quantities (essentially K0,K,Kr,y,e) of the system at a set of temperatures {Tn}, at various length scales (l, no greater than a length scale corresponding to the choice of lmax), and at some pressure (corresponding to the choice of Cc, the critical dimensionless vortex-loop-core energy).
   * For a set of temperatures {Tn}, the program starts from some initial state (l=0,K0=K0(Cc,Tn),K=K0,Kr=K0,y=y(K0,Cc),e=0) (depending on each temperature Tn and the choice of Cc) at a smallest length scale (l=0 and a=a0), it integrates to larger length scales (incrementing by Dl using the Runge-Kutta technique) up to some limiting length scale (see while loop, where l<lmax and K<5), and it estimates and records the state at that final length scale (l=lf,K0=K0(Cc,Tn),K=Kf,Kr=Kr(K,l),y=yf,e=ef).
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
const double B     = 4.0*PICU/3.0;
const double THETA = 0.6;

// Program parameters/inputs //
const int    N        = 3;         //  Selecting quantities to be calculated (2 -> K and y/G; 3 -> K, y/G, and e)
const int    Op       = 1;         //  Option (1 -> hand-picked temp's; else temperature-spread formulas)
const double lmax     = 100.0;     //  1  //  maximum of the length scale l (potentially approximate)
const int    lsteps   = 10000;     //  1250, 5000, 100000 //  lsteps = lmax/dl = 10/0.0001 = 100,000  (make sure this is divisible by ldatamax)
const int    lpts     = lsteps+1;  //  from l=0 to l=lmax, inclusive (incl. the boundary-condition-enforcing pnt, there are lpts+1 "l" pnts)
const int    ldatamax = 250;       //  50,100 //  max number of data points recorded per temperature examined
const double a0       = 1.0;       //  a0 in units of a0 (!)
const double a03      = a0*a0*a0;  //  a0 to the third power
const double a06      = a03*a03;   //  a0 to the sixth power
const double Kthresh  = 5.0;       //  ...



//  Function Prototypes  //
//=========================================================================//

// rk4((*f)(),*x,y[],h,n);
// EqRecRel(x,z,dzdx,n);
void rk4(void (*f)(double, double*, double*, unsigned int), double *x, double y[], double h, unsigned n);
void EqRecRel(double x, double z[], double dzdx[], unsigned n);



//  Function Definitions  //
//=========================================================================//

int main(){
	// Main function definitions //
	int    s,i,j, notdone=true, TempLabelOp;
	double dblsteps=lsteps, dl=lmax/dblsteps;
	double x, l[lpts];
	double Tfrac, tempv;
	double z[N], K[lpts], G[lpts], e[lpts];
	double K0, G0, y0;
	double Kr[lpts], Cc, K0c;
	FILE   *outfile;
	char   *filename;

	// Progress through sets of values of Cc and K0c, where a data file will be produced for each set //
	for(s=1; s<=1; s++){
		switch(s){
			// from plots_earlyequil/plots09/fit_Pdep_Cc_results.dat and 3Dvlt_K0cFind.c //
			case 1:	Cc = 1.10500797392886;	K0c = 0.295614012972922;	break;  // P = 0.050;
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
			case 22:	Cc = 0.1;	K0c = 1.52472890244101;		break;*/
			/* from vlt_K0cFind.c using vltRecRel with B' = 2B: // redo them */
			/* from vlt_K0cFind.c using vltRecRel with B' = B/2: // redo them */
		}
		printf(         "\n%s%g\t%s%g%s\n\n", " ---  Cc = ",Cc,"K0c = ",K0c,"  ---");

		// Select temperatures (calculations for each temp are independent), where a data file will be produced for each temp //
		i=0;
		while(notdone==true){
			switch(Op){
				case 1:	// A few hand-picked temperatures //
					switch(i){
						case 0:  Tfrac=1.00;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 1:  Tfrac=0.90;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 2:  Tfrac=0.80;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 3:  Tfrac=0.70;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 4:  Tfrac=0.60;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 5:  Tfrac=0.50;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 6:  Tfrac=0.40;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 7:  Tfrac=0.30;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 8:  Tfrac=0.20;  tempv=1.0-Tfrac;  i++;            TempLabelOp=0;  break;
						case 9:  Tfrac=0.10;  tempv=1.0-Tfrac;  notdone=false;  TempLabelOp=0;  break;
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
			if(TempLabelOp==0)  asprintf(&filename, "3Dvlt_constT_%g_Cc%3.2f_lmax%g_dl%g_Op%i.dat",  Tfrac,Cc,lmax,dl,Op);
			if(TempLabelOp==1)  asprintf(&filename, "3Dvlt_constTv_%g_Cc%3.2f_lmax%g_dl%g_Op%i.dat", tempv,Cc,lmax,dl,Op);
			outfile = fopen(filename,"w");  //  E.g., "3Dvlt_constT.out"
			printf(     "\n\n# Filename: %s\n", filename);
			fprintf(outfile,"# Filename: %s\n", filename);
			fprintf(outfile,"# Source: 3Dvlt_constT.c\n");
			fprintf(outfile,"# Source version: %s\n", "8 (2011 Jun 30 - ...)");
			fprintf(outfile,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g, B=%21.21g, THETA=%g\n", PI,PISQ,PICU,B,THETA);
			fprintf(outfile,"# Parameter values: N=%i, Tfrac=%g, tempv=%g, Cc=%g, K0c=%g, lmax=%g, lsteps=%i, (dl=%g), lpts=%i, ldatamax=%i, a0=%g, a03=%g, a06=%g\n", N,Tfrac,tempv,Cc,K0c,lmax,lsteps,dl,lpts,ldatamax,a0,a03,a06);
			fprintf(outfile,"# Option values: Op=%i\n", Op);
			fprintf(outfile,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","Kr/K0","G","e");
			printf(         "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","Kr/K0","G","e");

			// Initialize quantities at smallest length-scale (l=0) //
			x = l[0] = 0.0;
			z[0]  = K[0] = K0 = K0c/(1.0-tempv);
			z[1]  = y0 = exp(-PISQ*K0*Cc);
			z[2]  = e[0] = 0.0;
			Kr[0] = K[0]*exp(-l[0]);
			G[0]  = z[1]*exp(-6.0*l[0])/a06;

			// Calculate state-variables at progressively larger length scales, up to coherence length, via Runge-Kutta method (rk4) //
			j=0;
			while(l[j]<lmax && K[j]<Kthresh){
				// Step out in length scale l, calculating next l, K, y, and e (and G and Kr) //
				rk4(EqRecRel, &x, z, dl, N);  // equil: rk4 3D (K,y,e), EqRecRel
				j++;
				l[j]  = l[j-1] + dl;
				K[j]  = z[0];
				G[j]  = z[1]*exp(-6.0*l[j])/a06;
				e[j]  = z[2];
				Kr[j] = K[j]*exp(-l[j]);
				// Print at most [ldatamax] (e.g., 50) data points to screen and output file //
				if(j%(lsteps/ldatamax)==0){
					printf(         "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", Tfrac, tempv, l[j], K[j], K[j]/K0, Kr[j]/K0, G[j], e[j]);
					fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", Tfrac, tempv, l[j], K[j], K[j]/K0, Kr[j]/K0, G[j], e[j]);
				}
			}

			fclose(outfile);
		}
	}
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



// Equilibrium (K,y,e) recursion relations //
void EqRecRel(double x, double z[], double dzdx[], unsigned n){
	dzdx[0] = z[0]-B*z[0]*z[0]*z[1];
	if( log(z[0]) > 0.0 ) // if ac' > a, set ac' = a (replace log(K[i]) with zero)
		dzdx[1] = z[1]*( 6.0 - PISQ*z[0] );
	else
		dzdx[1] = z[1]*( 6.0 - PISQ*z[0]*(1.0-THETA*log(z[0])) );
	dzdx[2] = -PI*z[1]*exp(-3.0*x);
}



//  Program Notes  //
//=========================================================================//
/*

...

*/
