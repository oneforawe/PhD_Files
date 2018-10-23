//  File Comments  //
//=========================================================================//

/* FILENAME: 3Dvlt_KstarFind.c
   * Taken from  http://www.gnu.org/software/gsl/manual/html_node/Root-Finding-Examples.html
   * See  demo_RootFinder3.c  */
/* COMPILE NOTES:
   * To compile, type:
     g++ -L/usr/local/include/gsl-1.15/roots/.libs 3Dvlt_KstarFind.c -lgsl -lgslcblas -lgslroots -lm
   * To run, type
     ./a.out
   * These command formats were found at
   *      http://www.gnu.org/software/gsl/manual/html_node/Compiling-and-Linking.html
   * and  http://www.gnu.org/software/gsl/manual/html_node/Linking-with-an-alternative-BLAS-library.html
   * etc
*/



//  Function Preparation and Definition  //
//=========================================================================//

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>


const double PI    = 3.14159265358979323846;
const double PISQ  = 9.86960440108935799230;
const double PICU  = 31.00627668029981620634;
const double A3    = 4.0*PICU/3.0;
const double THETA = 0.6;

// Latest results from this program //
///////////////////////////////////////////////////////////////////////////////
const double KstarA = 4.146766416946759293;  //  4.146766416946759293
const double ystarA = 1.0/(A3*KstarA);       //  0.0058331356032128735
const double KstarB = 0.387508189712343409;  //  0.387508189712343409
const double ystarB = 1.0/(A3*KstarB);       //  0.0624210054576019
const double KstarC = 0.0;
const double ystarC = 0.0;
const double KstarD = 6.0/PISQ;        //                  0.607927101854026652
const double ystarD = 1.0/(A3*KstarD); // = 1.0/(8.0*PI) = 0.039788735772973836
const double nu     = -2.0  /  ( 1.0 - sqrt(1.0+4.0*PISQ*KstarB*(1.0-THETA*(1.0+log(KstarB)))) );  // = 0.671688352279845358
const double alpha  = 2.0 - 3.0*nu;  // = -0.015065056839536073
///////////////////////////////////////////////////////////////////////////////

struct CriticalFlowEqn_params{
	double constant, theta;
};

double CriticalFlowEqn (double x, void *params);


int
main (void){
	int status;
	int iter = 0, max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
//	double r = 4.0, r_expected = 4.15;   //  for Kstar "A"
//	double x_lo = 4.0, x_hi = 4.3;       //  for Kstar "A"
	double r = 0.0, r_expected = 0.388;  //  for Kstar "B"
	double x_lo = 0.3, x_hi = 0.4;       //  for Kstar "B"
	gsl_function F;
	struct CriticalFlowEqn_params params = {6.0, 0.6};

	F.function = &CriticalFlowEqn;
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
		status = gsl_root_test_interval (x_lo, x_hi, 0, 1.0e-6);  //  for Kstar "2"
//		status = gsl_root_test_interval (x_lo, x_hi, 0, 1.0e-6);  //  for Kstar "1"

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

	printf ("The converged `brent' root above gives the fixed-point value, Kstar: %.18f\n", r);
	printf ("For the record, here are all four of the various fixed-point value-sets (and critical exponents), previously calculated:\n");
	printf ("\tFixed point 'A':\n");
	printf ("\t\tKstarA = %.18f\n", KstarA );
	printf ("\t\tystarA = 1/(A3*KstarA) = %.18f\n", 1.0/(A3*KstarA) );
	printf ("\tFixed point 'B':\n");
	printf ("\t\tKstarB = %.18f\n", KstarB );
	printf ("\t\tystarB = 1/(A3*KstarB) = %.18f\n", 1.0/(A3*KstarB) );
	printf ("\tFixed point 'C' is trivial:\n");
	printf ("\t\tKstarC = 0\n");
	printf ("\t\tystarC = 0\n");
	printf ("\tFixed point 'D' can easily be calculated:\n");
	printf ("\t\tKstarD = 6/PISQ = %.18f\n", 6.0/PISQ );
	printf ("\t\tystarD = 1/(A3*KstarD) = 1/(8*PI) = %.18f\n", 1.0/(8.0*PI) );
	printf ("\tCritical exponents 'B':\n");
	printf ("\t\tnu    =  -2  /  ( 1 - sqrt(1+4*PISQ*KstarB*(1-THETA*(1+log(KstarB)))) )  =  %.18f\n", nu);
	// nu =  -2.0  /  ( 1.0 - sqrt(1.0+4.0*x) )
	// x  = PISQ*KstarB*( 1.0-THETA*(1.0+log(KstarB)) )
	printf ("\t\talpha =  2 - 3*nu  =  %.18f\n", alpha);

	return status;
}



double
CriticalFlowEqn (double x, void *params){
	struct CriticalFlowEqn_params *p 
	= (struct CriticalFlowEqn_params *) params;

	double constant  = p->constant;
	double theta     = p->theta;

	return constant - PISQ*x*(1.0-theta*log(x));
}
