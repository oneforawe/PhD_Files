/*  File Comments  */
/*=========================================================================*/

/* FILENAME: 2DKT_FPquench_adapt_gsl.c */
/* VERSION: 1 */
/* AUTHOR: ? (found at http://www.gnu.org/software/gsl/manual/gsl-ref.html ) */
/* DESCRIPTION:
   * Runge-Kutta integration of differential equations using adaptive step size.
*/
/* EXT FILES: none */
/* COMPILE NOTES:
   * Based on suggestions here ( http://www.gnu.org/software/gsl/manual/html_node/Compiling-and-Linking.html )
     and here ( http://www.gnu.org/software/gsl/manual/html_node/Shared-Libraries.html )
     follow similar instructions, such as the following (which works for me):
   * To compile, type (without the quotes) "g++ -Wall -I/usr/include -c test_gsl_rungekutta_clear.c"
   * To link to the gsl libraries, type    "g++ -L/usr/lib -lm -lgsl -lgslcblas test_gsl_rungekutta_clear.o"
   * To run, type                          "./a.out".
   * NOTE: Be sure...
*/



/*  Function Preparation  */
/*=========================================================================*/
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

//#define dim 2
const int dim = 2;



/*  Function prototypes  */
/*=========================================================================*/
//  int func(t,y[],f[],*params);
//  int jac(t,y[],*dfdy,dfdt[],*params);
int func(double t, const double y[], double f[], void *params);
int jac(double t, const double y[], double *dfdy, double dfdt[], void *params);



/*  Function definitions  */
/*=========================================================================*/

int main(void){
	FILE *outfile;
	const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45; // gsl_odeiv_step_rk8pd;
	double h = 1.8e-4;

	gsl_odeiv_step    * s = gsl_odeiv_step_alloc (T, dim);
	gsl_odeiv_control * c = gsl_odeiv_control_y_new (h, 0.0);
	  // control: maintain an absolute accuracy (error per step) of h and relative error per step of 0.0 in the function values y.
	gsl_odeiv_evolve  * e = gsl_odeiv_evolve_alloc (dim);

	double mu = 10;
	gsl_odeiv_system sys = {func, jac, dim, &mu};

	double t = 0.0, t1 = 100.0;
	//double h = 1e-6;
	double y[dim] = {2.5, 1.0};

	/* Prepare output file, print identification and values */
	sprintf(filename, "2DKT_FPquench_adapt_gsl_lmax%i_dl%g_tmax%g_dt0%g.dat", lmax, dl, tmax, dt0);
	outfile = fopen(filename,"w");  //  E.g., "2DKT_FPquench_adapt_gsl.out" or "2DKT_FPquench_adapt_gsl.dat"
	fprintf(outfile,"# Filename: %s\n",filename);
	fprintf(outfile,"# Source: 2DKT_FPquench_adapt_gsl.c\n");
	fprintf(outfile,"# t (time)  x = y[0]    x' = y[1]\n");
	printf(         "# t (time)  x = y[0]    x' = y[1]\n");

	while(t<t1){
		fprintf(outfile,"%.5e %.5e %.5e\n", t, y[0], y[1]);
		printf(         "%.5e %.5e %.5e\n", t, y[0], y[1]);
		int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, t1, &h, y);
		if(status!=GSL_SUCCESS)  break;
	}

	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);
	return 0;
}



int func(double t, const double y[], double f[], void *params){
	double mu = *(double *)params;
	f[0] = y[1];
	f[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1);
	// (two equations => system of dimension dim = 2)
	return GSL_SUCCESS;
}



int jac(double t, const double y[], double *dfdy, double dfdt[], void *params){
	double mu = *(double *)params;
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, dim, dim);
	gsl_matrix * m = &dfdy_mat.matrix;
	gsl_matrix_set (m, 0, 0, 0.0);
	gsl_matrix_set (m, 0, 1, 1.0);
	gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1]-1.0);
	gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0]-1.0));
	/* Reset? */
	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	return GSL_SUCCESS;
}
