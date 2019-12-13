//  File Comments  //
//=========================================================================//

/* FILENAME: 3Dvlt_macro.c */
/* VERSION: 8 (2011 Jun 30 - 2011 Aug 17)
     This program was formerly called vlt_ThermStates.c. */
/* AUTHOR: Andrew Forrester <aforrester@ucla.edu> */
/* DESCRIPTION:
   * This C program calculates macroscopic (or system-scale) thermodynamic properties of a superliquid helium film at various temperatures, given a certain chosen system size lmax, and puts the results into one output file.
   * This program calculates macroscopic equilibrium thermodynamic (or quasi-thermo-static) properties of a 2D system of superfluid (a thin layer of liquid 4He), using the Kosterlitz-Thouless (KT) theory, which could be considered a "vortex pair theory" (vpt, as opposed to the 3D vortex loop theory, vlt), hence the name "2Dvpt_macro.c".  By macroscopic, we mean that the quantities calculated are only recorded at the largest length scale, lmax (the size of the system), rather than recording all the renormalized values from the smallest scale (l=0 or a=a0) all the way up to lmax.
   * This program starts from a particular initial condition at the smallest length scale a0 and estimates the conditions at larger length scales for particular temperatures using differential recursion relations (obtained from the KT vortex-pair He-transition theory) and a 4th-order Runge-Kutta method.
   * The system is a thin film of liquid helium (4He) at atmospheric pressure, under or near the critical transition temperature Tc = Tkt separating normal fluid and superfluid phases.
   * The properties calculated are the quantities K, K/K0=sigma_s/sigma_0, y, and Gamma (not e?) as functions of temperature T (or Tfrac or tempv), which depend on the maximum length scale lmax (or amax) chosen.  See the Program Notes below the code for explanation of constants, variables, and functions.
*/
/* EXT FILES: none */
/* COMPILE NOTES:
   * To compile, type "g++ -lm 3Dvlt_macro.c" without the quotes; then to run, type "./a.out".
   * NOTE: Be sure lsteps is divisible by ldatamax.
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
const int    Op       = 3;         //  Option (1 -> hand-picked temp's; else temperature-spread formulas)
const double lmax     = 3.0;     //  1  //  maximum of the length scale l (potentially approximate)
const int    lsteps   = 3000;    //  1250, 5000, 100000 //  lsteps = lmax/dl = 10/0.0001 = 100,000  (make sure this is divisible by ldatamax)
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
	double K0, y0;
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

		switch(Op){  //  TempLabelOp -> 0 or 1 for labeling the file with "T" or "Tv", respectively
			case 1:  TempLabelOp=0;  break;
			case 2:  TempLabelOp=1;  break;
			case 3:  TempLabelOp=0;  break;
			case 4:  TempLabelOp=0;  break;
			case 5:  TempLabelOp=1;  break;
		} // TempLabelOp, and the file label (T or Tv) is meant to indicate which quantity is best used for plotting the data T/Tc = T/Tkt or Tv = tempv = 1-T/Tc = 1-T/Tkt

		// Prepare output file, print identification and values //
		if(TempLabelOp==0)  asprintf(&filename, "3Dvlt_macro_T_Cc%3.2f_lmax%g_dl%g_Op%i.dat",  Cc,lmax,dl,Op);
		if(TempLabelOp==1)  asprintf(&filename, "3Dvlt_macro_Tv_Cc%3.2f_lmax%g_dl%g_Op%i.dat", Cc,lmax,dl,Op);
		outfile = fopen(filename,"w");  //  E.g., "3Dvlt_constT.out"
		printf(     "\n\n# Filename: %s\n", filename);
		fprintf(outfile,"# Filename: %s\n", filename);
		fprintf(outfile,"# Source: 3Dvlt_macro.c\n");
		fprintf(outfile,"# Source version: %s\n", "8 (2011 Jun 30 - ...)");
		fprintf(outfile,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g, B=%21.21g, THETA=%g\n", PI,PISQ,PICU,B,THETA);
		fprintf(outfile,"# Parameter values: N=%i, Cc=%g, K0c=%g, lmax=%g, lsteps=%i, (dl=%g), lpts=%i, ldatamax=%i, a0=%g, a03=%g, a06=%g\n", N,Cc,K0c,lmax,lsteps,dl,lpts,ldatamax,a0,a03,a06);
		fprintf(outfile,"# Option values: Op=%i\n", Op);
		fprintf(outfile,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","Kr/K0","G","e");
		printf(         "\n%s%g\t%s%g%s\n\n", " ---  Cc = ",Cc,"K0c = ",K0c,"  ---");
		printf(         "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "T/Tc","1-T/Tc","l","K","K/K0","Kr/K0","G","e");

		/* Perform calculations for various temperatures (calculations for each temp are independent) */
		i=0;
		while(notdone==true){
			switch(Op){
				case 1:	// A few hand-picked temperatures //
					switch(i){
						case 0:   Tfrac=1.00;  tempv=1.0-Tfrac;  i++;            break;
						case 1:   Tfrac=0.90;  tempv=1.0-Tfrac;  i++;            break;
						case 2:   Tfrac=0.80;  tempv=1.0-Tfrac;  i++;            break;
						case 3:   Tfrac=0.70;  tempv=1.0-Tfrac;  i++;            break;
						case 4:   Tfrac=0.60;  tempv=1.0-Tfrac;  i++;            break;
						case 5:   Tfrac=0.50;  tempv=1.0-Tfrac;  i++;            break;
						case 6:   Tfrac=0.40;  tempv=1.0-Tfrac;  i++;            break;
						case 7:   Tfrac=0.30;  tempv=1.0-Tfrac;  i++;            break;
						case 8:   Tfrac=0.20;  tempv=1.0-Tfrac;  i++;            break;
						case 9:   Tfrac=0.10;  tempv=1.0-Tfrac;  notdone=false;  break;
					}
					break;
				case 2:	// Decrement below Tc (staying very close to Tc) //
					// from T/Tc = 1-1e-8 to T/Tc = 0.992
					tempv = 0.0+pow(10,-8.0+0.1*i);  Tfrac = 1.0-tempv;
					if(i==59) notdone=false;
					i++;  break;
				case 3:	// Increment from near zero to above Tc //
					// from T/Tc = 0.206 to T/Tc = 0.9990;  then from T/Tc = 1.0012 to T/Tc = 2
					if(i<30){ tempv = 0.0+pow(10,-0.1-0.1*i);      Tfrac = 1.0-tempv; }
					else    { tempv = 0.0-pow(10,-5.0+0.1*(i-9));  Tfrac = 1.0-tempv; }
					if(i==59) notdone=false;
					i++;  break;
				case 4:	// Increment from near zero to just below Tc //
					// from T/Tc = 0.206 to T/Tc = 0.9991
					tempv = 0.0+pow(10,-0.1-0.05*i);  Tfrac = 1.0-tempv;
					if(i==59) notdone=false;
					i++;  break;
				case 5:	// Decrement below Tc (staying extremely close to Tc) //
					// T/Tc = 1-1e-20  to  T/Tc ~ 1-1e-16
					tempv = 0.0+pow(10,-20.0+0.1*i);  Tfrac = 1.0-tempv;
					// earlier: T/Tc = 1-1e-12  to  T/Tc ~ 1-1e-6
					if(i==59) notdone=false;
					i++;  break;
			}

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
				// Bug Check: Print at most [ldatamax] (e.g., 50) data points to screen and output file //
				//if(j%(lsteps/ldatamax)==0){
				//	printf(         "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", Tfrac, tempv, l[j], K[j], K[j]/K0, Kr[j]/K0, G[j], e[j]);
				//	fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", Tfrac, tempv, l[j], K[j], K[j]/K0, Kr[j]/K0, G[j], e[j]);
				//}
			}

			/* Record data values */
			printf(         "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", Tfrac, tempv, l[j], K[j], K[j]/K0, Kr[j]/K0, G[j], e[j]);
			fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", Tfrac, tempv, l[j], K[j], K[j]/K0, Kr[j]/K0, G[j], e[j]);
			fflush(outfile);
		}

		fclose(outfile);
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
