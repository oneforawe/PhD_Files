/* filename: ring2.c */
/* description: This program calculates thermodynamic properties of a system of liquid helium at different temperatures and length scales.  It starts from a particular initial condition and estimates the conditions at larger length scales and lower temperatures using differential recursion relations (obtained from a vortex-loop He-transition theory) and a 4th-order Runge-Kutta method.
       The system is (an infinite bulk of? a spherical volume of?) liquid helium (4He) at atmospheric pressure, slightly under the critical transition temperature Tc separating normal fluid and superfluid.
       The properties calculated are the quantities K, y, and F (and Kr/K0 and a0/Kr) as functions of temperature T (or tempv) and length scale l (or a).  See below the code for explanation of constants, variables, and functions.
       The filename ring refers to the vortex rings (or loops) in this theory.
       (The program ring2.c is my second working nontrivial variation on the ring.c program that Dr. Gary Williams wrote.)
 */
/* ext files: derk2.c */
/* author: Andrew Forrester <aforrester@ucla.edu>  */

#include <stdio.h>
#include <math.h>

#define IN
#define OUT
#define INOUT
#define PI    3.14159265
#define PISQ  9.8696044
#define a0    2.53
#define K0c   0.3091469984216758
#define C     1.03
#define A     4.0*PI*PI*PI/3.0
#define N     3



/*  Function prototypes  */
/*=========================================================================*/
//   derk(func(),calculated,l,Dl,n)
//   ringrecrel(l,z)

extern void derk(void (*func)(double l, double z[], double dzdl[]), double calculated[], double *lPtr, double Dl, unsigned n);

void ringrecrel(double l, double z[], double dzdl[]);



/*  Function definitions  */
/*=========================================================================*/

int main(){
	double K0,Kr,a;
	int i,j=0;
	int Dexp=100;  // 100 corresponds to D = 6.8e33 meters (huge system)
	double tempv, l, estimated[N];
	FILE *outfile;

	/* Print headings for data to screen and output file */
	outfile = fopen("ring2_B_Dexp100.dat","w");  //  "ring2.out"
	fprintf(outfile,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
		"T/Tc","1-T/Tc","l","a","K","y","F","Kr","Kr/K0","a0/Kr");
	printf(         "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
		"T/Tc","1-T/Tc","l","a","K","y","F","Kr","Kr/K0","a0/Kr");

	/* Loop through different temperatures (calculations for each temp are independent) */
	for(i=0; i<60; i++){
		/* (Option A) Decrement temperature value (further below Tc)
		tempv = 0.0+pow(10,-8.0+0.1*i); */
		/* (Option B) Increment temperature value (near zero to above Tc) */
		if(i<30) tempv = 0.0+pow(10,-0.1-0.1*i);
		else     tempv = 0.0-pow(10,-5.0+0.1*(i-9));

		/* Initialize temperature-dependent values */
		K0 = K0c/(1.0-tempv);
		estimated[0] = K0;
		estimated[1] = exp(-PISQ*K0*C);

		/* Initialize lengthscale and energy values */
		l = 0;
		estimated[2] = 0;

		/* Progress to larger lengthscales via Runge-Kutta method */
		while(l<Dexp && estimated[0]<5.0){
			/* Calculate next estimated state, incrementing lengthscale */
			derk(IN ringrecrel, INOUT estimated, INOUT &l, IN 0.0001, IN N);  // Dl=0.0001, n=3
			/* Print out intermediate numbers
			a = a0*exp(l);
			Kr = z[0]*a0/a;
			j++;
				fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
					1-tempv, tempv, l, a, estimated[0], estimated[1], estimated[2], Kr, Kr/K0, a0/Kr);  // writes to a temporary file (not yet to outfile, see fflush below)
				printf(         "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
					1-tempv, tempv, l, a, estimated[0], estimated[1], estimated[2], Kr, Kr/K0, a0/Kr);
			}
			if(j>100000)  j=0; */
		}

		/* Update dependent values, if not using "Print out intermediate numbers" block of code */
		a = a0*exp(l);
		Kr = estimated[0]*a0/a;

		/* Record data values */
		fprintf(outfile,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
			1-tempv, tempv, l, a, estimated[0], estimated[1], estimated[2], Kr, Kr/K0, a0/Kr);  // writes to a temporary file (not yet to outfile, see fflush below)
		fflush(outfile);  // forces write-to-file (outfile) from the temporary file (see fprintf above), so if an error interrupts the execution of the code at least some data will be recorded (for bug tracking)
		printf(         "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
			1-tempv, tempv, l, a, estimated[0], estimated[1], estimated[2], Kr, Kr/K0, a0/Kr);
	}

	fclose(outfile);
	return 0;
}



void ringrecrel(double l, double z[], double dzdl[]){
	dzdl[0] = z[0]-A*z[1]*z[0]*z[0];
	dzdl[1] = z[1]*(6.0-PISQ*z[0]*(1.0-0.6*log(z[0])));
	dzdl[2] = -PI*z[1]*exp(-3.0*l);
}
