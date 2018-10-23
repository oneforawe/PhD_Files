//  File Comments  //
//=========================================================================//

/* FILENAME: 2Dvpt_K0cFind.c */
/* VERSION: 2 (2012 Jul 31 - ...)
            Based off of 3Dvlt_KstarFind.c
            First version, combining code from 3Dvlt_K0cFind.c and 2Dvpt_macro.c */
/* FILENAME: 3Dvlt_KstarFind.c
   * Taken from  http://www.gnu.org/software/gsl/manual/html_node/Root-Finding-Examples.html
   * See  demo_RootFinder3.c  */
/* COMPILE NOTES:
   * To compile, type:
     g++ -L/usr/local/include/gsl-1.15/roots/.libs 2Dvpt_K0cFind.c -lgsl -lgslcblas -lgslroots -lm
   * To run, type
     ./a.out
   * These command formats were found at
   *      http://www.gnu.org/software/gsl/manual/html_node/Compiling-and-Linking.html
   * and  http://www.gnu.org/software/gsl/manual/html_node/Linking-with-an-alternative-BLAS-library.html
   * etc
*/

/* LATEST RESULTS:
   * K0c = 0.747852421772060461
     y0c = 0.000622972891299253652
*/



//  Function Preparation and Definition  //
//=========================================================================//

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>


const double PI       = 3.14159265358979323846;
const double PISQ     = 9.86960440108935799230;
const double PICU     = 31.00627668029981620634;
const double A2       = 4.0*PICU;
const double Kstar    = 2.0/PI; // = 0.636619772367581382
const double Vmfactor = 1.0;    //  factor to adjust the Villain model
// Villain model:  U  =  Vmfactor*PISQ*K0*kB*T + (kappa^2*sigma^r_s/(2*PI))*ln(a/a0)
//                 U0 =  Vmfactor*PISQ*K0*kB*T


struct CriticalFlowIntersectEqn_params {
	double xf;
};

double CriticalFlowIntersectEqn (double x, void *params);


int
main (void){
	int status;
	int iter = 0, max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r = 0.0, r_expected = 0.748; // can change these
	double x_lo = 0.7, x_hi = 0.8;      // can change these
	gsl_function F;
	struct CriticalFlowIntersectEqn_params params = {0.636619772367581382};

	F.function = &CriticalFlowIntersectEqn;
	F.params = &params;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc (T);
	gsl_root_fsolver_set (s, &F, x_lo, x_hi);

	printf("\nUsing %s method:\n", gsl_root_fsolver_name (s));

	printf ("%5s [%20s, %20s] %20s %21s %20s\n",
		"iter", "lower", "upper", "root", 
		"err", "err(est)");

	do{
		iter++;
		status = gsl_root_fsolver_iterate (s);
		r = gsl_root_fsolver_root (s);
		x_lo = gsl_root_fsolver_x_lower (s);
		x_hi = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (x_lo, x_hi, 0, 1.0e-6);  //  what're these again?

		if (status == GSL_SUCCESS)
		printf ("Converged:\n");

		printf ("%5d [%.18f, %.18f] %.18f %+.18f %.18f\n",
			iter, x_lo, x_hi,
			r, r - r_expected,
			x_hi - x_lo);
	}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free (s);
	printf ("\n");

	printf ("The converged `brent' root above gives the bare critical value, K0c: %.18f\n", r);
	printf ("\tThus, y0c = exp(-PISQ*K0c) = %.18g\n", exp(-PISQ*r) );
	printf ("For the record, here is the special, critical fixed-point value-set (within the line of fixed points along y=0):\n");
	printf ("\t(Critical) fixed point:\n");
	printf ("\t\tKstar = 2/PI = %.18f\n", 2.0/PI );
	printf ("\t\tystar = 0\n");

	return status;
}



double
CriticalFlowIntersectEqn (double x, void *params){
	struct CriticalFlowIntersectEqn_params *p 
	= (struct CriticalFlowIntersectEqn_params *) params;

	double xf  = p->xf;

	return exp(-Vmfactor*PISQ*x) + (4.0/A2)*(1.0/xf - 1.0/x) + (2.0*PI/A2)*log(xf/x);
}



// equilibrium (K,y) recursion relations //
void EqRecRel(double x, double z[2], double dzdx[2], unsigned n){
	dzdx[0] = -A2*z[0]*z[0]*z[1];
	dzdx[1] = z[1]*(4.0-2.0*PI*z[0]);
}



//  Program Notes  //
//=========================================================================//
/*

See further documentation (UCLAresearchNotes.pdf, UCLAresearchPrograms.pdf, vlt_ThermStates.c, ring1.c) for elaboration.

Dl is the increment in l.
n is the number of thermodynamic variables (n=3 given thermv[i], i=1,2,3).

DthrmvDl[i].arr[0] = dK/dl
DthrmvDl[i].arr[1] = dy/dl
DthrmvDl[i].arr[2] = de/dl

*/
