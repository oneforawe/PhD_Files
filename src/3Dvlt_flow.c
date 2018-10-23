//  File Comments  //
//=========================================================================//

/* FILENAME: 3Dvlt_flow.c */
/* VERSION: 1 (2012 Jul 19 - ...)
            Based on 3Dvlt_K0cFind3_Plot.c v2 (or v1 including the if-then statement in the recursion relations)
            Now using if-then statement in recursion relations to account for ac not surpassing a
            Now using a GSL (GNU Science Library) subroutine to find the intersection of the exponential and the trajectory.
            Now using a method that is more informed by the autonomous differential equations analysis, using dy/dK. */
/* AUTHOR: Andrew Forrester <aforrester@ucla.edu> */
/* DESCRIPTION:
   * ... (look in 3Dvlt_K0cFind.c to mimick notes.)

   Inputs:
   * A file called "3Dvlt_K0cFind3_Inputn.dat" where "n" matches the number written below.  E.g., vlt_CcInput3.dat.
     In the file there should be a list of values for Cc.  Cc can be any positive number, presumably.
   * ... (look in 3Dvlt_K0cFind.c to mimick notes.)

   Output:
   * ... (look in 3Dvlt_K0cFind.c to mimick notes.)
*/
/* EXT FILES: vlt_CcInputn.dat */
/* COMPILE NOTES:
   * To compile, first prepare a new input file, then change the input/output filenames below,
     and then type "g++ -lm 3Dvlt_flow.c" without the quotes; then to run, type "./a.out".
*/
/* PROGRAM IDEAS:
   * When writing a program like this, one should not just calculate the values of interest but also the estimated error in each step summed up to give a total error.
   * We may need to use an arbitrary-precision package to make sure we're getting all the digits we need/want.
   * To make this code more efficient, realize that there is really only one trajectory that determines all of the K0c's.  The list of Cc's can be loaded all at once, arranged into two lists (one for the positive-K-direction trajectory portion, and the other for the negative-K-direction trajectory portion), and each K0c can be calculated along the way as the trajectory is traversed.
   * One could also make the code more efficient by turning some of the blocks of code into subroutine functions, called several times.
*/

/* RESULTS:
*/



//  Function Preparation  //
//=========================================================================//

// Standard routine header files //
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Constants and labels //
const double PI     = 3.14159265358979323846;
const double PISQ   = 9.86960440108935861883;
const double PICU   = 31.00627668029981620634;

// Major option combinations
// Op  A    B      C        D         E
//     123  12345  1(2345)  1234(56)  1234  //  MAJOR OPTION 1: Any flow type, any starting point (including fixed pts B,D & eigenvector paths)
                                            //                  (If using itlog flow, start from a physically relevant point to get meaningful data)
//     2    3      23       7         1234  //  MAJOR OPTION 2: Finding K0c's (w/log flow, starting at fixed pt B, using eigenvector paths)

// Program parameters/inputs //
const int    OpA  = 1;   //  Option A   Set recursion relation behavior
                         //              1 -> "itlog flow": use if-then statement to limit ac at a (one fixed point, B) -- "if-then log flow" of EqRecRel1;
                         //              2 -> "w/log flow": no if-then, using logarithm so ac may exceed a (two fixed points, A,B) -- "with log flow" of EqRecRel2;
                         //              3 -> "nolog flow": no if-then, no logarithm so ac = a the whole time (two fixed pts, C,D) -- "no log flow" of EqRecRel3.
const int    OpB  = 1;   //  Option B   Set starting point in a flow trajectory
                         //              1 -> arbitrary pt: start at some user-defined point hand-entered below;
                         //              2 -> fixed pt A:   start at w/log-flow fixed-point "1", the lower-right spiral fixed point;
                         //              3 -> fixed pt B:   start at w/log-flow fixed-point "2", the upper-left saddle fixed point;
                         //              4 -> fixed pt C:   start at nolog-flow fixed-point "1", the origin unstable node fixed point;
                         //              5 -> fixed pt D:   start at nolog-flow fixed-point "2", the upper-right saddle fixed point.
const int    OpC  = 1;   //  Option C   Set/Determine initial trajectory path/axis (not yet determining forward or backward direction, ds or dl)
                         //              1 -> if using any starting pt, path determined by flow trajectory at starting pt (cutoff arranged in code below);
                         //              2 -> if using fixed pt B, traverse the attractive path (lambda<0) to/from fixed pt;
                         //              3 -> if using fixed pt B, traverse the repulsive path  (lambda>0) to/from fixed pt;
                         //              4 -> if using fixed pt D, traverse the attractive path (lambda<0) to/from fixed pt;
                         //              5 -> if using fixed pt D, traverse the repulsive path  (lambda>0) to/from fixed pt.
const int    OpD  = 4;   //  Option D   Set first step direction (forward or backward along flow, or -/+ in K-direction)
                         //              1 -> first step direction use fwd flow trajectory at initial pt (const dl);
                         //              2 -> first step direction use bwd flow trajectory at initial pt (const dl);
                         //              3 -> first step direction use fwd flow trajectory at initial pt (limited ds);
                         //              4 -> first step direction use bwd flow trajectory at initial pt (limited ds);
                         //              5 -> first step direction from fixed point (const ds) leftward  (- in K-direction);
                         //              6 -> first step direction from fixed point (const ds) rightward (+ in K-direction);
                         //              7 -> first step direction from fixed point (const ds) determined from Cc and destination curve (y0 = exp(-PISQ*Cc*K0)) to find K0c.
const int    OpE  = 4;   //  Option E   Set (2nd and later steps) direction and const dl steps or limited ds steps (lengths in l-space or K-y-space)
                         //              1 -> second+ steps use fwd traj const dl step size;
                         //              2 -> second+ steps use rev traj const dl step size;
                         //              3 -> second+ steps use fwd traj limited ds step size;
                         //              4 -> second+ steps use rev traj limited ds step size.
const double A3     = 4.0*PICU/3.0;
const double THETA  = 0.6;

const double KstarA = 4.146766416946759293;  //  4.146766416946759293
const double ystarA = 1.0/(A3*KstarA);       //  0.0058331356032128735
const double KstarB = 0.387508189712343409;  //  0.387508189712343409
const double ystarB = 1.0/(A3*KstarB);       //  0.0624210054576019
const double KstarC = 0.0;
const double ystarC = 0.0;
const double KstarD = 6.0/PISQ;        //                  0.607927101854026652
const double ystarD = 1.0/(A3*KstarD); // = 1.0/(8.0*PI) = 0.039788735772973836

const int    N      = 2;
const double initialdl       = 1.0e-5;  // first step size, if using flow
const double trajectorydl    = 1.0e-5;  // step size (if using dl) taken in trajectory, using scaling/recursion equations - (second and later steps if using eigenvector direction for first step)
const double initialds       = 1.0e-5;  // first step size, if using eigenvector direction
const double trajectorydsMax = 1.0e-5;  // step size (if using ds) taken in trajectory, using scaling/recursion equations - (second and later steps)
const double uncertainty     = 1.0e-16;  // loop, calculating K0c 2+ times, decreasing dl until change in K0c is within uncertainty
const int    dataskipdl = floor(1.0/trajectorydl/100.0);     // The first 10 points in the trajectory are recorded and then (given const dl)   every "dataskipdl"th point is recorded (for plotting)
const int    dataskipds = floor(1.0/trajectorydsMax/100.0);  // The first 10 points in the trajectory are recorded and then (given limited ds) every "dataskipds"th point is recorded (for plotting)
const int    Max_steps  = 1000000;
	// Max_steps =  1000000 (for const dl, left-downward) or 4000000;
	//              1500000 (for const dl, left-upward);
	//              2000000 (for const dl, right-downward);
	//             11000000 (for const dl, right-upward);
	//              4000000 (for something...)
const double FudgeFactor = 1.0e0;



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
	int     i; // Cc-case counter (when taking Cc values from input file and finding K0c's)
	double  Kstar = -10.0;
	double  ystar = -10.0;
	double  z[N], K, y;
	double  oldK, oldy;
	double  dK, dy, deltaK;
	double  ds, dl, dlnew;
	double  Deltas=0.0, Deltal=0.0; // path-length-s (in K-y-space) and path-length-l
	int     dataskipdlAdj, dataskipdsAdj;
	double  l=0.0; // nothing happens to l anyway since it's irrelevant (Deltal is relevant, though)
	int     initialKdir=0, flowdir=0;
	int     nextK0cfound;
	int64_t count_steps=0;
	int     count_passes=0;
	char    trash[4];
	double  Cc, K0c=0.0, oldK0c=-100.0;
	double  c=0.0, Z=0.0, lambda=0.0;
	FILE    *infile,*outfile1,*outfile2;
	char    *outfilename1,*outfilename2;
	void    (*EqRecRel)(double x, double z[], double dzdx[], unsigned n);

	switch(OpA){ // Set recursion relation behavior
		case 1:  // "itlog flow": use if-then statement to limit ac at a (one fixed point, B) -- "if-then log flow"
			 EqRecRel = EqRecRel1;
			 break;
		case 2:  // "w/log flow": no if-then, using logarithm so ac may exceed a (two fixed points, A,B) -- "with log flow"
			 EqRecRel = EqRecRel2;
			 break;
		case 3:  // "nolog flow": no if-then, no logarithm so ac = a the whole time (two fixed pts, C,D) -- "no log flow"
			 EqRecRel = EqRecRel3;
			 break;
	}

	switch(OpB){ // Set starting point in a flow trajectory
		case 1:  // arbitrary pt: start at some user-defined point hand-entered below
			 // NOTE: If using itlog flow (OpA=1), start from a physically relevant point to get meaningful data !!!

			 //Cc = 0.40;  // Max_steps=1000000 (dl=1.0e-5, lmax<=10)
			 //K0c = 0.57534342964893814;
			 Cc = 1.10;   // Max_steps=1000000 (dl=1.0e-5, lmax<=10)
			 K0c = 0.29646893900678273;

			 K = K0c/2.2; //= K0 @ T=2.2*Tc
			 //K = K0c/2.1; //= K0 @ T=2.1*Tc
			 //K = K0c/2.0; //= K0 @ T=2.0*Tc
			 //K = K0c/1.9; //= K0 @ T=1.9*Tc
			 //K = K0c/1.8; //= K0 @ T=1.8*Tc
			 //K = K0c/1.7; //= K0 @ T=1.7*Tc
			 //K = K0c/1.6; //= K0 @ T=1.6*Tc
			 //K = K0c/1.5; //= K0 @ T=1.5*Tc
			 //K = K0c/1.4; //= K0 @ T=1.4*Tc
			 //K = K0c/1.3; //= K0 @ T=1.3*Tc
			 //K = K0c/1.2; //= K0 @ T=1.2*Tc
			 //K = K0c/1.1; //= K0 @ T=1.1*Tc
			 //K = K0c/1.05;//= K0 @ T=1.05*Tc
			 //K = K0c; //= K0 @ T=Tc
			 //K = K0c/0.95;//= K0 @ T=0.95*Tc
			 //K = K0c/0.9; //= K0 @ T=0.9*Tc
			 //K = K0c/0.8; //= K0 @ T=0.8*Tc
			 //K = K0c/0.7; //= K0 @ T=0.7*Tc
			 //K = K0c/0.6; //= K0 @ T=0.6*Tc
			 //K = K0c/0.5; //= K0 @ T=0.5*Tc
			 //K = K0c/0.4; //= K0 @ T=0.4*Tc
			 //K = K0c/0.3; //= K0 @ T=0.3*Tc
			 //K = K0c/0.2; //= K0 @ T=0.2*Tc
			 //K = K0c/0.1; //= K0 @ T=0.1*Tc
			 y = exp(-PISQ*K*Cc); //= y0

			 //K = KstarA;
			 //y = ystarA;  //1000000 and 1000000
			 //K = 2.1;
			 //y = 1.0/(A3*K);  //1000000 and 1000000
			 //K = 0.92226247542;
			 //y = 1.0/(A3*K);  //1000000 and 1000000
			 //K = KstarB;
			 //y = ystarB;  //1000000 and 1000000
			 //K = 1.0/(A3*0.15);
			 //y = 0.15;  //1000000 and 1000000

			 //K = KstarD;
			 //y = 0.2;  //1000000 and 1000000
			 //K = KstarD;
			 //y = 0.125;  //1000000 and 1000000
			 //K = KstarD;
			 //y = 1.5*ystarD;  //1000000 and 1000000
			 //K = KstarD;
			 //y = 0.5*ystarD;  //1000000 and 1000000
			 //K = KstarD;
			 //y = 0.125*ystarD;  //1000000 and 1000000

			 // From using K0cFind3:
			 //K = 0.387508189712343409409;
			 //y = -0.075;  //1000000 and 1000000
			 //K = 4.1449909483947875;
			 //y = 0.25709876948489824-0.1; //100000 and 100000
			 //K = 4.1449909483947875;
			 //y = 0.25709876948489824-0.2; //1000000 and 100000
			 //K = 32.6;
			 //y = 0.39; //1000000 
			 //K = 32.6;
			 //y = 0.29; //1000000 
			 //K = 32.6;
			 //y = 0.19; //1000000 
			 //K = 32.6;
			 //y = 0.09; //1000000 
			 //K = 32.6;
			 //y = -0.01; //1000000 and 1000
			 //K = 32.6;
			 //y = -0.11; //1000000 and 1000
			 //K = 32.6;
			 //y = -0.21; //1000000 and 1000
			 //K = 32.6;
			 //y = -0.31; //1000000 and 1000
			 //K = 32.6;
			 //y = -0.41; //1000000 and 1000
			 //K = 32.6;
			 //y = -0.51; //1000000 and 1000
			 //K = 32.6;
			 //y = -0.61; //1000000 and 1000
			 //K = 32.6;
			 //y = -0.71; //1000000 and 1000
			 //K = 0.38750;
			 //y = 0.06238; //
			 //K = 0.38750;
			 //y = 0.06246; //
			 //K = 0.387375;
			 //y = 0.06244; //
			 //K = 0.38765;
			 //y = 0.06240; //
			 //K = 0.38750;
			 //y = 0.06244; //100000 and 100000
			 //K = 0.387575;
			 //y = 0.06242; //100000 and 150000
			 //K = 0.38755;
			 //y = 0.06242; //100000 and 100000
			 //K = 0.38745;
			 //y = 0.06242; //100000 and 150000
			 //K = 0.38750;
			 //y = 0.06240; //100000 and 100000
			 //K = 0.4;
			 //y = 0.025; //150000 and 1200000
			 //K = 0.4;
			 //y = 0.09; //500000 and 1000000
			 //K = 0.15;
			 //y = 0.06; //1000000 and 300000
			 //K = 1.0e-4;
			 //y = 1.0e-75; //1500000
			 //K = 1.0;
			 //y = 1.0e-20; //1500000 and 1000000
			 //K = 100.0;
			 //y = 1.0e-40; //100000
			 //K = 40.0;
			 //y = 1.0e-40; //1000000
			 //K = 32.0;
			 //y = 1.0e-40; //1000000
			 //K = 25.0;
			 //y = 1.0e-40; //1000000
			 //K = 1.0e-5;
			 //y = 1.0e-20; //800000
			 //K = 1.0e-13;
			 //y = 1.0e-10; //800000
			 //K = 1.0e-10;
			 //y = 1.0e13; // bah! couldn't get any good data for a flow curve
			 //K = 1.0;
			 //y = 3.0; //500000
			 //K = 1.0;
			 //y = 1.0;
			 //K = 0.07;
			 //y = -0.001;
			 //K = 0.04;
			 //y = -0.001;
			 //K = 0.05;
			 //y = -0.001;
			 //K = 0.1;
			 //y = -0.001;
			 //K = 0.01;
			 //y = -0.001;
			 //K = 0.005;
			 //y = -0.001;
			 //K = 3.0;
			 //y = 0.05;
			 //K = 0.25;
			 //y = 0.0125;
			 //K = 0.25;
			 //y = 0.25;
			 //K = 7.5;
			 //y = 0.002;
			 break;
		case 2:  // fixed pt A:   start at w/log-flow fixed-point "1", the lower-right spiral fixed point
			 K = Kstar = KstarA; //                   4.146766416946759293
			 y = ystar = ystarA; // 1.0/(A3*KstarA) = 0.0058331356032128735
			 // ? c = ystar * PISQ * ( 1.0-THETA*(1.0+log(Kstar)) );
			 // ? Z = Kstar * PISQ * ( 1.0-THETA*(1.0+log(Kstar)) );
			 break;
		case 3:  // fixed pt B:   start at w/log-flow fixed-point "2", the upper-left saddle fixed point
			 K = Kstar = KstarB; //                   0.387508189712343409
			 y = ystar = ystarB; // 1.0/(A3*KstarB) = 0.0624210054576019
			 c = ystar * PISQ * ( 1.0-THETA*(1.0+log(Kstar)) );
			 Z = Kstar * PISQ * ( 1.0-THETA*(1.0+log(Kstar)) );
			 break;
		case 4:  // fixed pt C:   start at nolog-flow fixed-point "1", the origin unstable node fixed point
			 K = Kstar = KstarC; // 0.0
			 y = ystar = ystarC; // 0.0
			 break;
		case 5:  // fixed pt D:   start at nolog-flow fixed-point "2", the upper-right saddle fixed point
			 K = Kstar = KstarD; // 6.0/PISQ =                       0.607927101854026652
			 y = ystar = ystarD; // 1.0/(A3*KstarD) = 1.0/(8.0*PI) = 0.039788735772973836
			 c = ystar * PISQ;
			 break;
	}
	z[0] = K;
	z[1] = y;

	switch(OpC){ // Set/Determine initial trajectory path/axis (not yet determining forward or backward direction, ds or dl)
		case 1:  // if using any starting pt, path determined by flow trajectory at starting pt (cutoff arranged in code below)
			 break;
		case 2:  // if using fixed pt B, traverse the attractive path (lambda<0) to/from fixed pt
			 lambda = 0.5*( -1.0 - sqrt(1.0+4.0*Z) );
			 break;
		case 3:  // if using fixed pt B, traverse the repulsive path  (lambda>0) to/from fixed pt
			 lambda = 0.5*( -1.0 + sqrt(1.0+4.0*Z) );
			 break;
		case 4:  // if using fixed pt D, traverse the attractive path (lambda<0) to/from fixed pt
			 lambda = -3.0;
			 break;
		case 5:  // if using fixed pt D, traverse the repulsive path  (lambda>0) to/from fixed pt
			 lambda = 2.0;
			 break;
	}

	switch(OpD){ // Set first step direction (forward or backward along flow, or -/+ in K-direction)
		case 1:  // first step direction use fwd flow trajectory at initial pt (const dl)
			 flowdir = +1;
			 break;
		case 2:  // first step direction use bwd flow trajectory at initial pt (const dl)
			 flowdir = -1;
			 break;
		case 3:  // first step direction use fwd flow trajectory at initial pt (limited ds)
			 flowdir = +1;
			 break;
		case 4:  // first step direction use bwd flow trajectory at initial pt (limited ds)
			 flowdir = -1;
			 break;
		case 5:  // first step direction from fixed point (const ds) leftward  (- in K-direction)
			 initialKdir = -1;
			 break;
		case 6:  // first step direction from fixed point (const ds) rightward (+ in K-direction)
			 initialKdir = +1;
			 break;
		case 7:  // first step direction from fixed point (const ds) determined from Cc and destination curve (y0 = exp(-PISQ*Cc*K0)) to find K0c
			 break; // determination made below
	}



	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//                                                                                                        //
	//  MAJOR OPTION 1: Any flow type, any starting point (including fixed pts B,D & eigenvector paths)       //
	//                  (If using itlog flow, start from a physically relevant point to get meaningful data)  //
	//                                                                                                        //
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if( OpD!=7 ){
		// Op  A    B      C        D         E
		//     123  12345  1(2345)  1234(56)  1234
		// Any flow type, any starting point (including fixed pts B,D & eigenvector paths)
		// (first step direction determined by options)
		// (cutoff arranged in code below with Max_steps)

		// Get ready to record trajectory points for plotting //
		asprintf(&outfilename2, "3Dvlt_flow_Op%i%i%i%i%i_K%4.3f_y%4.3f_step%i.dat", OpA,OpB,OpC,OpD,OpE,K,y,Max_steps);
//		asprintf(&outfilename2, "3Dvlt_flow_Op%i%i%i%i%i_K%6.5f_y%6.5f_step%i.dat", OpA,OpB,OpC,OpD,OpE,K,y,Max_steps);
		outfile2 = fopen(outfilename2,"w");  //  E.g., "3Dvlt_flow.out"
		fprintf(outfile2,"# Filename: %s\n", outfilename2);
		fprintf(outfile2,"# Source: 3Dvlt_flow.c\n");
		fprintf(outfile2,"# Source version: %s\n", "1 (2012 Jul 19 - ...)");
		fprintf(outfile2,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g\n", PI,PISQ,PICU);
		fprintf(outfile2,"# Parameter values: OpA=%i, OpB=%i, OpC=%i, OpD=%i, OpE=%i, A3=%21.21g, THETA=%g, N=%i, initialdl=%g, trajectorydl=%g, initialds=%g, trajectorydsMax=%g, uncertainty=%g, dataskipdl=%i, dataskipds=%i\n", OpA, OpB, OpC, OpD, OpE, A3, THETA, N, initialdl, trajectorydl, initialds, trajectorydsMax, uncertainty, dataskipdl, dataskipds);
		fprintf(outfile2,"# Option-dependent values: K=%21.21g, y=%21.21g, Kstar=%21.21g, ystar=%21.21g, c=%21.21g, Z=%21.21g, lambda=%21.21g, initialKdir=%i, flowdir=%i, Max_steps=%i, FudgeFactor=%g\n", K, y, Kstar, ystar, c, Z, lambda, initialKdir, flowdir, Max_steps, FudgeFactor);
		fprintf(outfile2,"# Initialized values: Deltas=%g, Deltal=%g, count_steps=%ld, count_passes=%i, K0c=%g, oldK0c=%g\n", Deltas, Deltal, count_steps, count_passes, K0c, oldK0c);
		fprintf(outfile2,"# Input-file parameter: Cc=%g\n", Cc);
		fprintf(outfile2,"# %s\t%s\t%s\t%s\t%s\n","step","K","y","path-length-l","path-length-s (in K-y-space)");
		printf(          "# %s\t%s\t%s\t%s\t%s\n","step","K","y","path-length-l","path-length-s (in K-y-space)");

		// Print initialized position for trajectory //
		printf(          "%ld\t%17.17g\t%17.17g\t%g\t%g\n",count_steps,K,y,Deltal,Deltas);
		fprintf(outfile2,"%ld\t%17.17g\t%17.17g\t%g\t%g\n",count_steps,K,y,Deltal,Deltas);
		fflush(outfile2);

		// Take first step //
		if(OpD==1||OpD==2){ // const dl
			oldK = K;
			oldy = y;
			dl = flowdir*initialdl;
			rk4(EqRecRel, &l, z, dl, N);  // equil: rk4 3D (K,y), EqRecRel
			K = z[0];
			y = z[1];
			count_steps++;
			Deltal += dl;
			Deltas += sqrt((K-oldK)*(K-oldK)+(y-oldy)*(y-oldy));
		}
		if(OpD==3||OpD==4){ // limited ds
			oldK = K;
			oldy = y;
			dl = flowdir*initialdl;
			dlnew = dl;
			rk4(EqRecRel, &l, z, dl, N);  // equil: rk4 3D (K,y), EqRecRel
			ds = sqrt((z[0]-oldK)*(z[0]-oldK)+(z[1]-oldy)*(z[1]-oldy));
			if(ds>trajectorydsMax){
				//printf("ds = %g  >  trajectorydsMax = %g\n", ds,trajectorydsMax);
				z[0] = oldK;
				z[1] = oldy;
				//dlnew = (trajectorydsMax/ds)*dl;
				dlnew = FudgeFactor*(trajectorydsMax/ds)*dl;
				rk4(EqRecRel, &l, z, dlnew, N);
				ds = sqrt((z[0]-oldK)*(z[0]-oldK)+(z[1]-oldy)*(z[1]-oldy));
				while(ds>trajectorydsMax){
					//printf("ds = %g  >  trajectorydsMax = %g\n", ds,trajectorydsMax);
					z[0] = oldK;
					z[1] = oldy;
					dlnew = dlnew/2.0;
					rk4(EqRecRel, &l, z, dlnew, N);
					ds = sqrt((z[0]-oldK)*(z[0]-oldK)+(z[1]-oldy)*(z[1]-oldy));
				}
			}
			K = z[0];
			y = z[1];
			count_steps++;
			Deltal += dlnew;
			Deltas += sqrt((K-oldK)*(K-oldK)+(y-oldy)*(y-oldy));
		}
		if(OpD==5||OpD==6){ // const ds (-/+ in K-direction)
			oldK = K;
			oldy = y;
			dK = initialKdir*(1.0/sqrt(1.0+c*c/lambda/lambda))*initialds;
			dy = -(c/lambda)*dK;
			K += dK;
			y += dy;
			z[0] = K;
			z[1] = y;
			count_steps++;
			// Deltal = can't update because the change in l is/may-be infinite, from the fixed point
			Deltas += initialds;
			//Deltas += sqrt((K-oldK)*(K-oldK)+(y-oldy)*(y-oldy));
		}
		// Print first step //
		fprintf(outfile2,"%ld\t%17.17g\t%17.17g\t%g\t%g\n",count_steps,K,y,Deltal,Deltas);
		printf(          "%ld\t%17.17g\t%17.17g\t%g\t%g\n",count_steps,K,y,Deltal,Deltas);
		fflush(outfile2);

		// Take additional steps //
		if(OpE==1||OpE==3) flowdir = +1;
		if(OpE==2||OpE==4) flowdir = -1;
		dl = flowdir*trajectorydl;
		if(OpE==1||OpE==2){ // const dl
			while(count_steps<Max_steps){ // cutoff arranged in code here !!
				oldK = K;
				oldy = y;
				rk4(EqRecRel, &l, z, dl, N);  // equil: rk4 3D (K,y), EqRecRel
				K = z[0];
				y = z[1];
				count_steps++;
				Deltal += dl;
				Deltas += sqrt((K-oldK)*(K-oldK)+(y-oldy)*(y-oldy));
				if( count_steps<11 || count_steps%dataskipdl==0 ){
					fprintf(outfile2,"%ld\t%17.17g\t%17.17g\t%g\t%g\n",count_steps,K,y,Deltal,Deltas);
					printf(          "%ld\t%17.17g\t%17.17g\t%g\t%g\n",count_steps,K,y,Deltal,Deltas);
					fflush(outfile2);
				}
			}
		}
		if(OpE==3||OpE==4){ // limited ds
			while(count_steps<Max_steps){ // cutoff arranged in code here !!
				oldK = K;
				oldy = y;
				dlnew = dl;
				rk4(EqRecRel, &l, z, dl, N);  // equil: rk4 3D (K,y), EqRecRel
				ds = sqrt((z[0]-oldK)*(z[0]-oldK)+(z[1]-oldy)*(z[1]-oldy));
				if(ds>trajectorydsMax){
					//printf("ds = %g  >  trajectorydsMax = %g\n", ds,trajectorydsMax);
					z[0] = oldK;
					z[1] = oldy;
					//dlnew = (trajectorydsMax/ds)*dl;
					dlnew = FudgeFactor*(trajectorydsMax/ds)*dl;
					rk4(EqRecRel, &l, z, dlnew, N);
					ds = sqrt((z[0]-oldK)*(z[0]-oldK)+(z[1]-oldy)*(z[1]-oldy));
					while(ds>trajectorydsMax){
						//printf("ds = %g  >  trajectorydsMax = %g\n", ds,trajectorydsMax);
						z[0] = oldK;
						z[1] = oldy;
						dlnew = dlnew/2.0;
						rk4(EqRecRel, &l, z, dlnew, N);
						ds = sqrt((z[0]-oldK)*(z[0]-oldK)+(z[1]-oldy)*(z[1]-oldy));
					}
				}
				K = z[0];
				y = z[1];
				count_steps++;
				Deltal += dlnew;
				Deltas += sqrt((K-oldK)*(K-oldK)+(y-oldy)*(y-oldy));
				if( count_steps<11 || count_steps%dataskipds==0 ){
					fprintf(outfile2,"%ld\t%17.17g\t%17.17g\t%g\t%g\n",count_steps,K,y,Deltal,Deltas);
					printf(          "%ld\t%17.17g\t%17.17g\t%g\t%g\n",count_steps,K,y,Deltal,Deltas);
					fflush(outfile2);
				}
			}
		}

		// (Bug checking) When checking approach to fixed point B //
		//printf("%ld\t%17.17g\t%17.17g\t%g\t%g\n",count_steps,KstarB,ystarB,Deltal,Deltas);

		// Close plotting output file //
		fclose(outfile2);
	}



	///////////////////////////////////////////////////////////////////////////////////////////////////
	//                                                                                               //
	//  MAJOR OPTION 2: Finding K0c's (w/log flow, starting at fixed pt B, using eigenvector paths)  //
	//                                                                                               //
	///////////////////////////////////////////////////////////////////////////////////////////////////
	if( OpD==7 ){
		// Op  A    B      C        D         E
		//     2    3      23       7         1234
		// Need to read in Cc's from infile and will spit out K0c's into outfile1,
		// along with spitting out points to plot in outfile2.
		// (first step direction determined by infile Cc)
		// (cutoff also determined by infile Cc)

		// Open/name input and output files, and print headings for data to screen and output file //
		infile  = fopen(        "3Dvlt_flow_Input1.dat","r");
		asprintf(&outfilename1, "3Dvlt_flow_Output1_OpE%i_uncert%g.dat", OpE,uncertainty);
		outfile1 = fopen(outfilename1,"w");  //  E.g., "3Dvlt_flow.out"
		fprintf(outfile1,"# Filename: %s\n", outfilename1);
		fprintf(outfile1,"# Source: 3Dvlt_flow.c\n");
		fprintf(outfile1,"# Source version: %s\n", "1 (2012 Jul 19 - ...)");
		fprintf(outfile1,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g\n", PI,PISQ,PICU);
		fprintf(outfile1,"# Parameter values: OpA=%i, OpB=%i, OpC=%i, OpD=%i, OpE=%i, A3=%21.21g, THETA=%g, N=%i, initialdl=%g, trajectorydl=%g, initialds=%g, trajectorydsMax=%g, uncertainty=%g, dataskipdl=%i, dataskipds=%i\n", OpA, OpB, OpC, OpD, OpE, A3, THETA, N, initialdl, trajectorydl, initialds, trajectorydsMax, uncertainty, dataskipdl, dataskipds);
		fprintf(outfile1,"# Option-dependent values: K=%21.21g, y=%21.21g, Kstar=%21.21g, ystar=%21.21g, c=%21.21g, Z=%21.21g, lambda=%21.21g, initialKdir=%i, flowdir=%i, Max_steps=%i\n", K, y, Kstar, ystar, c, Z, lambda, initialKdir, flowdir, Max_steps);
		fprintf(outfile1,"# Initialized values: Deltas=%g, Deltal=%g, count_steps=%ld, count_passes=%i, K0c=%g, oldK0c=%g\n", Deltas, Deltal, count_steps, count_passes, K0c, oldK0c);
		fprintf(outfile1,"# %s\t%s\n","Cc","K0c");
		printf(            "%s\t%s\n","Cc","K0c");
		fflush(outfile1);

		// Skip first line in input file to ready data retrieval //
		fscanf(infile,"%s",trash);  // Pass "Cc"

		// Find K0c for each (Cc,K01) pair in the input file //
		i = 0;
		while(!feof(infile)){

			// Acquire value for Cc from input file //
			fscanf(infile,"%lf",&Cc);
			if(feof(infile)) break;
			// Increment Cc,K0c data-set/case count (for output file) //
			i++;

			// Until K0c is found within the desired precision, ... //
			count_passes = 0;
			while( fabs(K0c-oldK0c) > uncertainty ){

				// Re-initialize starting point (for second and later passes of while loops) //
				count_passes++;
				count_steps = 0;
				K = Kstar;
				y = ystar;
				Deltal = 0.0;
				Deltas = 0.0;

				// Get ready to (re)record trajectory points for plotting (throwing away faulty, imprecise passes) //
				asprintf(&outfilename2, "3Dvlt_flow_Op%i%i%i%i%i_Cc%3.2f_uncert%g.dat", OpA,OpB,OpC,OpD,OpE,Cc,uncertainty);
				outfile2 = fopen(outfilename2,"w");  //  E.g., "3Dvlt_flow.out"
				fprintf(outfile2,"# Filename: %s\n", outfilename2);
				fprintf(outfile2,"# Source: 3Dvlt_flow.c\n");
				fprintf(outfile2,"# Source version: %s\n", "1 (2012 Jul 19 - ...)");
				fprintf(outfile2,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g\n", PI,PISQ,PICU);
				fprintf(outfile2,"# Parameter values: OpA=%i, OpB=%i, OpC=%i, OpD=%i, OpE=%i, A3=%21.21g, THETA=%g, N=%i, initialdl=%g, trajectorydl=%g, initialds=%g, trajectorydsMax=%g, uncertainty=%g, dataskipdl=%i, dataskipds=%i\n", OpA, OpB, OpC, OpD, OpE, A3, THETA, N, initialdl, trajectorydl, initialds, trajectorydsMax, uncertainty, dataskipdl, dataskipds);
				fprintf(outfile2,"# Option-dependent values: K=%21.21g, y=%21.21g, Kstar=%21.21g, ystar=%21.21g, c=%21.21g, Z=%21.21g, lambda=%21.21g, initialKdir=%i, flowdir=%i, Max_steps=%i\n", K, y, Kstar, ystar, c, Z, lambda, initialKdir, flowdir, Max_steps);
				fprintf(outfile2,"# Initialized values: Deltas=%g, Deltal=%g, count_steps=%ld, count_passes=%i, K0c=%g, oldK0c=%g\n", Deltas, Deltal, count_steps, count_passes, K0c, oldK0c);
				fprintf(outfile2,"# Input-file parameter: Cc=%g\n", Cc);
				fprintf(outfile2,"# Pass: count_passes=%i (to reach desired precision in K0c)\n", count_passes);
				fprintf(outfile2,"# %s\t%s\t%s\t%s\t%s\n","step","K","y","path-length-l","path-length-s (in K-y-space)");
				fflush(outfile2);

				// Print initialized position for trajectory //
				fprintf(outfile2,"%ld\t%17.17g\t%17.17g\t%g\t%g\n",count_steps,K,y,Deltal,Deltas);
				printf(          "%ld\t%17.17g\t%17.17g\t%g\t%g\n",count_steps,K,y,Deltal,Deltas);
				fflush(outfile2);

				// Find which direction (+/-) to travel in K-direction //
				if( y < exp(-PISQ*Cc*K) )
					initialKdir = +1;
				if( y > exp(-PISQ*Cc*K) )
					initialKdir = -1;
				if( y== exp(-PISQ*Cc*K) ){
					K0c = Kstar;
					fprintf(outfile1,"%17.17g\t%17.17g\n",Cc,K0c);
					printf(          "%17.17g\t%17.17g\n",Cc,K0c);
					fflush(outfile1);
					continue;
				}

				// Take first step in trajectory (if not already done finding K0c) //
				oldK = K;
				oldy = y;
				dK = initialKdir*(1.0/sqrt(1.0+c*c/lambda/lambda))*initialds;
				dy = -(c/lambda)*dK;
				K += dK;
				y += dy;
				z[0] = K;
				z[1] = y;
				count_steps++;
				// Deltal = can't update because the change in l is/may-be infinite, from the fixed point
				Deltas += initialds;
				//Deltas += sqrt((K-oldK)*(K-oldK)+(y-oldy)*(y-oldy));
				fprintf(outfile2,"%ld\t%17.17g\t%17.17g\t%g\t%g\n",count_steps,K,y,Deltal,Deltas);
				printf(          "%ld\t%17.17g\t%17.17g\t%g\t%g\n",count_steps,K,y,Deltal,Deltas);
				fflush(outfile2);

				// Check for K0c and take second+ steps in trajectory (if K0c not found) //
				if(OpE==1||OpE==2){ // const dl
					if(OpE==1) flowdir = +1; // this option is not being allowed at the moment since the code will cycle infinitely
					if(OpE==2) flowdir = -1;
					//dl = flowdir*trajectorydl;
					dl = flowdir*trajectorydl/pow(2.0,(count_passes-1));
					dataskipdlAdj = dataskipdl*pow(2.0,(count_passes-1));
					// dataskipdsAdj = dataskipds*pow(2.0,(count_passes-1)); don't need this here
//					if(count_passes==1){
//						dl = flowdir*trajectorydl;
//						dataskipdlAdj = dataskipdl;
//						dataskipdsAdj = dataskipds;
//					}
//					else{
//						dl = flowdir*trajectorydl*sqrt(uncertainty/fabs(K0c-oldK0c))/2.0;
//						dataskipdlAdj = dataskipdl*2.0/sqrt(uncertainty/fabs(K0c-oldK0c));
//						dataskipdsAdj = dataskipds*2.0/sqrt(uncertainty/fabs(K0c-oldK0c));
//					}

					// Take steps in trajectory, checking for interception (finding K0c) //
					nextK0cfound = false;
					if(initialKdir>0){
						while(nextK0cfound==false){
							if( y > exp(-PISQ*Cc*K) ){
								oldK0c = K0c;
								deltaK = (exp(-PISQ*Cc*oldK)-oldy)/((y-oldy)/(K-oldK)+PISQ*Cc*exp(-PISQ*Cc*oldK));
								K0c = oldK + deltaK;
								nextK0cfound = true;
								oldK = K;
								oldy = y;
								K = K0c;
								y = exp(-PISQ*Cc*K0c);
								count_steps++;
								// Deltal = should update somehow...
								Deltas += sqrt((K-oldK)*(K-oldK)+(y-oldy)*(y-oldy));
								fprintf(outfile2,"%ld\t%17.17g\t%17.17g\t%g\t%g\n",                  count_steps,K,y,Deltal,Deltas);
								printf(          "%ld\t%17.17g\t%17.17g\t%g\t%g\tpass=%i (Cc=%g) diff=%g goal=%g\n",count_steps,K,y,Deltal,Deltas,count_passes,Cc,fabs(K0c-oldK0c),uncertainty);
								fflush(outfile2);
							}
							else{
								oldK = K;
								oldy = y;
								rk4(EqRecRel, &l, z, dl, N);  // equil: rk4 3D (K,y), EqRecRel
								K = z[0];
								y = z[1];
								count_steps++;
								Deltal += dl;
								Deltas += sqrt((K-oldK)*(K-oldK)+(y-oldy)*(y-oldy));
								if( count_steps<11 || count_steps%dataskipdlAdj==0 ){
									fprintf(outfile2,"%ld\t%17.17g\t%17.17g\t%g\t%g\n",                  count_steps,K,y,Deltal,Deltas);
									printf(          "%ld\t%17.17g\t%17.17g\t%g\t%g\tpass=%i (Cc=%g) diff=%g goal=%g\n",count_steps,K,y,Deltal,Deltas,count_passes,Cc,fabs(K0c-oldK0c),uncertainty);
									fflush(outfile2);
								}
							}
						}
					}
					if(initialKdir<0){
						while(nextK0cfound==false){
							if( y < exp(-PISQ*Cc*K) ){
								oldK0c = K0c;
								deltaK = (exp(-PISQ*Cc*K)-y)/((y-oldy)/(K-oldK)+PISQ*Cc*exp(-PISQ*Cc*K));
								K0c = K + deltaK;
								nextK0cfound = true;
								oldK = K;
								oldy = y;
								K = K0c;
								y = exp(-PISQ*Cc*K0c);
								count_steps++;
								// Deltal = should update somehow...
								Deltas += sqrt((K-oldK)*(K-oldK)+(y-oldy)*(y-oldy));
								fprintf(outfile2,"%ld\t%17.17g\t%17.17g\t%g\t%g\n",                  count_steps,K,y,Deltal,Deltas);
								printf(          "%ld\t%17.17g\t%17.17g\t%g\t%g\tpass=%i (Cc=%g) diff=%g goal=%g\n",count_steps,K,y,Deltal,Deltas,count_passes,Cc,fabs(K0c-oldK0c),uncertainty);
								fflush(outfile2);
							}
							else{
								oldK = K;
								oldy = y;
								rk4(EqRecRel, &l, z, dl, N);  // equil: rk4 3D (K,y), EqRecRel
								K = z[0];
								y = z[1];
								count_steps++;
								Deltal += dl;
								Deltas += sqrt((K-oldK)*(K-oldK)+(y-oldy)*(y-oldy));
								if( count_steps<11 || count_steps%dataskipdlAdj==0 ){
									fprintf(outfile2,"%ld\t%17.17g\t%17.17g\t%g\t%g\n",                  count_steps,K,y,Deltal,Deltas);
									printf(          "%ld\t%17.17g\t%17.17g\t%g\t%g\tpass=%i (Cc=%g) diff=%g goal=%g\n",count_steps,K,y,Deltal,Deltas,count_passes,Cc,fabs(K0c-oldK0c),uncertainty);
									fflush(outfile2);
								}
							}
						}
					}
				}
				if(OpE==3||OpE==4){ // limited ds
					if(OpE==3) flowdir = +1; // this option is not being allowed at the moment since the code will cycle infinitely
					if(OpE==4) flowdir = -1;
					//dl = flowdir*trajectorydl;
					dl = flowdir*trajectorydl/pow(2.0,(count_passes-1));
					//dataskipdlAdj = dataskipdl*pow(2.0,(count_passes-1)); don't need this here
					dataskipdsAdj = dataskipds*pow(2.0,(count_passes-1));

					// Take steps in trajectory, checking for interception (finding K0c) //
					nextK0cfound = false;
					if(initialKdir>0){
						while(nextK0cfound==false){
							if( y > exp(-PISQ*Cc*K) ){
								oldK0c = K0c;
								deltaK = (exp(-PISQ*Cc*oldK)-oldy)/((y-oldy)/(K-oldK)+PISQ*Cc*exp(-PISQ*Cc*oldK));
								K0c = oldK + deltaK;
								nextK0cfound = true;
								oldK = K;
								oldy = y;
								K = K0c;
								y = exp(-PISQ*Cc*K0c);
								count_steps++;
								// Deltal = should update somehow...
								Deltas += sqrt((K-oldK)*(K-oldK)+(y-oldy)*(y-oldy));
								fprintf(outfile2,"%ld\t%17.17g\t%17.17g\t%g\t%g\n",                  count_steps,K,y,Deltal,Deltas);
								printf(          "%ld\t%17.17g\t%17.17g\t%g\t%g\tpass=%i (Cc=%g) diff=%g goal=%g\n",count_steps,K,y,Deltal,Deltas,count_passes,Cc,fabs(K0c-oldK0c),uncertainty);
								fflush(outfile2);
							}
							else{
								oldK = K;
								oldy = y;
								dlnew = dl;
								rk4(EqRecRel, &l, z, dl, N);  // equil: rk4 3D (K,y), EqRecRel
								ds = sqrt((z[0]-oldK)*(z[0]-oldK)+(z[1]-oldy)*(z[1]-oldy));
								if(ds>trajectorydsMax){
									//printf("ds = %g  >  trajectorydsMax = %g\n", ds,trajectorydsMax);
									z[0] = oldK;
									z[1] = oldy;
									//dlnew = (trajectorydsMax/ds)*dl;
									dlnew = FudgeFactor*(trajectorydsMax/ds)*dl;
									rk4(EqRecRel, &l, z, dlnew, N);
									ds = sqrt((z[0]-oldK)*(z[0]-oldK)+(z[1]-oldy)*(z[1]-oldy));
									while(ds>trajectorydsMax){
										//printf("ds = %g  >  trajectorydsMax = %g\n", ds,trajectorydsMax);
										z[0] = oldK;
										z[1] = oldy;
										dlnew = dlnew/2.0;
										rk4(EqRecRel, &l, z, dlnew, N);
										ds = sqrt((z[0]-oldK)*(z[0]-oldK)+(z[1]-oldy)*(z[1]-oldy));
									}
								}
								K = z[0];
								y = z[1];
								count_steps++;
								Deltal += dlnew;
								Deltas += ds;
								//Deltas += sqrt((K-oldK)*(K-oldK)+(y-oldy)*(y-oldy));
								if( count_steps<11 || count_steps%dataskipdsAdj==0 ){
									fprintf(outfile2,"%ld\t%17.17g\t%17.17g\t%g\t%g\n",                  count_steps,K,y,Deltal,Deltas);
//									printf(          "%ld\t%17.17g\t%17.17g\t%g\t%g\tpass=%i (Cc=%g) diff=%g goal=%g\n",count_steps,K,y,Deltal,Deltas,count_passes,Cc,fabs(K0c-oldK0c),uncertainty);
								}
								fflush(outfile2);
							}
						}
					}
					if(initialKdir<0){
						while(nextK0cfound==false){
							if( y < exp(-PISQ*Cc*K) ){
								oldK0c = K0c;
								deltaK = (exp(-PISQ*Cc*K)-y)/((y-oldy)/(K-oldK)+PISQ*Cc*exp(-PISQ*Cc*K));
								K0c = K + deltaK;
								nextK0cfound = true;
								oldK = K;
								oldy = y;
								K = K0c;
								y = exp(-PISQ*Cc*K0c);
								count_steps++;
								// Deltal = should update somehow...
								Deltas += sqrt((K-oldK)*(K-oldK)+(y-oldy)*(y-oldy));
								fprintf(outfile2,"%ld\t%17.17g\t%17.17g\t%g\t%g\n",                  count_steps,K,y,Deltal,Deltas);
								printf(          "%ld\t%17.17g\t%17.17g\t%g\t%g\tpass=%i (Cc=%g) diff=%g goal=%g\n",count_steps,K,y,Deltal,Deltas,count_passes,Cc,fabs(K0c-oldK0c),uncertainty);
								fflush(outfile2);
							}
							else{
								oldK = K;
								oldy = y;
								dlnew = dl;
								rk4(EqRecRel, &l, z, dl, N);  // equil: rk4 3D (K,y), EqRecRel
								ds = sqrt((z[0]-oldK)*(z[0]-oldK)+(z[1]-oldy)*(z[1]-oldy));
								if(ds>trajectorydsMax){
									//printf("ds = %g  >  trajectorydsMax = %g\n", ds,trajectorydsMax);
									z[0] = oldK;
									z[1] = oldy;
									//dlnew = (trajectorydsMax/ds)*dl;
									dlnew = FudgeFactor*(trajectorydsMax/ds)*dl;
									rk4(EqRecRel, &l, z, dlnew, N);
									ds = sqrt((z[0]-oldK)*(z[0]-oldK)+(z[1]-oldy)*(z[1]-oldy));
									while(ds>trajectorydsMax){
										//printf("ds = %g  >  trajectorydsMax = %g\n", ds,trajectorydsMax);
										z[0] = oldK;
										z[1] = oldy;
										dlnew = dlnew/2.0;
										rk4(EqRecRel, &l, z, dlnew, N);
										ds = sqrt((z[0]-oldK)*(z[0]-oldK)+(z[1]-oldy)*(z[1]-oldy));
									}
								}
								K = z[0];
								y = z[1];
								count_steps++;
								Deltal += dlnew;
								Deltas += ds;
								//Deltas += sqrt((K-oldK)*(K-oldK)+(y-oldy)*(y-oldy));
								if( count_steps<11 || count_steps%dataskipdsAdj==0 ){
									fprintf(outfile2,"%ld\t%17.17g\t%17.17g\t%g\t%g\n",                  count_steps,K,y,Deltal,Deltas);
									printf(          "%ld\t%17.17g\t%17.17g\t%g\t%g\tpass=%i (Cc=%g) diff=%g goal=%g\n",count_steps,K,y,Deltal,Deltas,count_passes,Cc,fabs(K0c-oldK0c),uncertainty);
									fflush(outfile2);
								}
							}
						}
					}
				}

				// Close plotting output file //
				fclose(outfile2);
			}

			// Print results to screen and output file //
			//fprintf(outfile1,"%17.17g\t%17.17g\n",Cc,K0c);
			//fprintf(outfile1,"case %i:\t\tCc = %1.15g;\tK0c = %1.15g;\tbreak;\n", i,Cc,K0c);
			fprintf(outfile1,"case %i:\t\tCc = %3.2g;\tK0c = %17.17g;\tbreak;\n", i,Cc,K0c);
			//fprintf(outfile1,"%17.17g\t%17.17g\t(steps=%ld, passes=%i, dK=%g)\n",Cc,K0c,count_steps,count_passes,dK);
			//fprintf(outfile1,"%17.17g\t%17.17g\t(steps=%ld, passes=%i, dK=%g)\n oldK0c = %17.17g\n    K0c = %17.17g\n   diff = %17.17g\n uncert = %17.17g\n\n",Cc,K0c,count_steps,count_passes,dK, oldK0c,K0c,K0c-oldK0c,uncertainty);
			printf(          "%17.17g\t%17.17g\t(steps=%ld, passes=%i)\n oldK0c = %17.17g\n    K0c = %17.17g\n   diff = %17.17g\n uncert = %17.17g\n\n",Cc,K0c,count_steps,count_passes, oldK0c,K0c,K0c-oldK0c,uncertainty);
			//printf(          "%17.17g\t%17.17g\t(steps=%ld, passes=%i, dK=%g)\n",Cc,K0c,count_steps,count_passes,dK); dK? What do I care about that?
			fflush(outfile1);

			// Get ready for next input and search //
			K0c = 0.0;
			oldK0c = -100.0;
			count_passes = 0;
		}

		// Close input and output files //
		fclose(infile);
		fclose(outfile1);
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
		//printf("K=%g\ty=%g\n",s[0],s[1]);
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

See further documentation (UCLAresearchNotes.pdf, UCLAresearchPrograms.pdf, vlt_ThermStates.c, ring1.c) for elaboration.

Earlier thought about Major Options:

// Major option combinations
// Op  A    B      C        D         E
//     1    1      1        1234      1234  //  MAJOR OPTION 1: arbitrary starting point, itlog flow (should start from a physically relevant point to be meaningful)
//     2    123    123      1234(567) 1234  //  MAJOR OPTION 2: w/log flow, any starting point (including fixed pt B & eigenvector paths, & finding K0c)
//     3    145    145      1234(56)  1234  //  MAJOR OPTION 3: nolog flow, any starting point (including fixed pt D & eigenvector paths)

*/
