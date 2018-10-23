//  File Comments  //
//=========================================================================//

/* FILENAME: 3Dvlt_flow.c */
/* VERSION: 1 (2012 Jul 31 - ...)
            Based on 3Dvlt_flow.c  */
/* AUTHOR: Andrew Forrester <aforrester@ucla.edu> */
/* DESCRIPTION:
   * ...

   Inputs:
   * ...

   Output:
   * ...
*/
/* EXT FILES: 2Dvpt_flow_Inputn.dat ?  */
/* COMPILE NOTES:
   * To compile, first prepare a new input file, then change the input/output filenames below,
     and then type "g++ -lm 2Dvpt_flow.c" without the quotes; then to run, type "./a.out".
*/
/* PROGRAM IDEAS:
   * When writing a program like this, one should not just calculate the values of interest but also the estimated error in each step summed up to give a total error.
   * We may need to use an arbitrary-precision package to make sure we're getting all the digits we need/want.
   * To make this code more efficient, realize that there is really only one trajectory that determines all of the K0c's.  The list of Cc's can be loaded all at once, arranged into two lists (one for the positive-K-direction trajectory portion, and the other for the negative-K-direction trajectory portion), and each K0c can be calculated along the way as the trajectory is traversed.
   * One could also make the code more efficient by turning some of the blocks of code into subroutine functions, called several times.
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
const int    OpA  = 4;   //  Option A   Set stepping direction (forward or reverse flow trajectory) and step-size type: const dl steps or limited ds steps (lengths in l-space or K-y-space)
                         //              1 -> fwd traj const dl step size;
                         //              2 -> rev traj const dl step size;
                         //              3 -> fwd traj limited ds step size;
                         //              4 -> rev traj limited ds step size.
const double A2     = 4.0*PICU;
const double Kstar  = 2.0/PI;  //  0.636619772368
const double ystar  = 0.0;     //  0.0
const double K0c    = 0.747852421772060461;
const int    N      = 2;
const double trajectorydl    = 1.0e-5;  // step size (if using dl) taken in trajectory, using scaling/recursion equations - (second and later steps if using eigenvector direction for first step)
const double trajectorydsMax = 1.0e-5;  // step size (if using ds) taken in trajectory, using scaling/recursion equations - (second and later steps)
const int    dataskipdl      = floor(1.0/trajectorydl/100.0);     // The first 10 points in the trajectory are recorded and then (given const dl)   every "dataskipdl"th point is recorded (for plotting)
const int    dataskipds      = floor(1.0/trajectorydsMax/100.0);  // The first 10 points in the trajectory are recorded and then (given limited ds) every "dataskipds"th point is recorded (for plotting)
const int    Max_steps       = 1000000;
const double FudgeFactor     = 1.0e0;
const double Vmfactor        = 1.0;    //  factor to adjust the Villain model
// Villain model:  U  =  factor*PISQ*K0*kB*T + (kappa^2*sigma^r_s/(2*PI))*ln(a/a0)
//                 U0 =  factor*PISQ*K0*kB*T



//  Function Prototypes  //
//=========================================================================//

// rk4((*f)(),*x,y[],h,n);
// EqRecRel(x,z,dzdx,n);
void rk4(void (*f)(double, double*, double*, unsigned int), double *x, double y[], double h, unsigned n);
void EqRecRel(double x, double z[], double dzdx[], unsigned n);



//  Function Definitions  //
//=========================================================================//

int main(void){
	double  z[N], K, y;
	double  oldK, oldy;
	double  dK, dy, deltaK;
	double  ds, dl, dlnew;
	double  Deltas=0.0, Deltal=0.0; // path-length-s (in K-y-space) and path-length-l
	double  l=0.0; // nothing happens to l anyway since it's irrelevant (Deltal is relevant, though)
	int     initialKdir=0, flowdir=0;
	int64_t count_steps=0;
	FILE    *outfile;
	char    *outfilename;

	// Set starting point in a flow trajectory //
	///////////////////////////////////////////////////////////////////////
	//K = Kstar;  // to create reference flow, from which to create the critical flow
	//y = 0.1;    // Max_steps=1000000
	//K = 10.622439796320707;        // critical flow to (Kstar,ystar)
	//y = -0.1+0.19496241950161122;  // OpA=3 Max_steps=20000000
	//K = 0.0034155620322175229;    // critical flow away from (Kstar,ystar)
	//y = -0.1+9.2270207538457303;  // OpA=4 Max_steps=20000000

	//K = K0c/1.075; //= K0 @ T=1.075*Tc  // Max_steps=5000000   1000000 (away)
	//K = K0c/1.050; //= K0 @ T=1.050*Tc  // Max_steps=5000000   1000000 (away)
	//K = K0c/1.025; //= K0 @ T=1.025*Tc  // Max_steps=5000000   1000000 (away)
	K = K0c;       //= K0 @ T=Tc        // Max_steps=10000000   1000000 (away)
	//K = K0c/0.975; //= K0 @ T=0.975*Tc  // Max_steps=5000000   1000000 (away)
	//K = K0c/0.950; //= K0 @ T=0.950*Tc  // Max_steps=5000000   1000000 (away)
	//K = K0c/0.925; //= K0 @ T=0.925*Tc  // Max_steps=5000000   1000000 (away)
	//K = K0c/0.900; //= K0 @ T=0.900*Tc  // Max_steps=1000000   1000000 (away)
	//K = K0c/0.850; //= K0 @ T=0.850*Tc  // Max_steps=1000000   1000000 (away)
	//K = K0c/0.800; //= K0 @ T=0.800*Tc  // Max_steps=1000000   1000000 (away)
	//K = K0c/0.750; //= K0 @ T=0.750*Tc  // Max_steps=1000000   1000000 (away)
	//K = K0c/0.700; //= K0 @ T=0.750*Tc  // Max_steps=1000000   1000000 (away)
	//K = K0c/0.650; //= K0 @ T=0.750*Tc  // Max_steps=1000000   1000000 (away)
	//K = K0c/0.600; //= K0 @ T=0.750*Tc  // Max_steps=1000000   1000000 (away)
	y = exp(-Vmfactor*PISQ*K); //= y0
	///////////////////////////////////////////////////////////////////////
	z[0] = K;
	z[1] = y;

	// Get ready to record trajectory points for plotting //
	asprintf(&outfilename, "2Dvpt_flow_Op%i_K%4.3f_y%4.3f_step%i.dat", OpA,K,y,Max_steps);
//	asprintf(&outfilename, "2Dvpt_flow_Op%i_K%6.5f_y%6.5f_step%i.dat", OpA,K,y,Max_steps);
	outfile = fopen(outfilename,"w");  //  E.g., "2Dvpt_flow.out"
	fprintf(outfile,"# Filename: %s\n", outfilename);
	fprintf(outfile,"# Source: 2Dvpt_flow.c\n");
	fprintf(outfile,"# Source version: %s\n", "1 (2012 Jul 31 - ...)");
	fprintf(outfile,"# Constant values: PI=%21.21g, PISQ=%21.21g, PICU=%21.21g\n", PI,PISQ,PICU);
	fprintf(outfile,"# Initialized values: Deltas=%g, Deltal=%g, count_steps=%ld\n", Deltas, Deltal, count_steps);
	fprintf(outfile,"# Option-dependent values: K=%21.21g, y=%21.21g, Kstar=%21.21g, ystar=%21.21g, flowdir=%i, Max_steps=%i, FudgeFactor=%g, Vmfactor=%g\n", K, y, Kstar, ystar, flowdir, Max_steps, FudgeFactor, Vmfactor);
	fprintf(outfile,"# Parameter values: OpA=%i, A2=%21.21g, N=%i, trajectorydl=%g, trajectorydsMax=%g, dataskipdl=%i, dataskipds=%i\n", OpA, A2, N, trajectorydl, trajectorydsMax, dataskipdl, dataskipds);
	fprintf(outfile,"# %s\t%s\t%s\t%s\t%s\n","step","K","y","path-length-l","path-length-s (in K-y-space)");
	printf(         "# %s\t%s\t%s\t%s\t%s\n","step","K","y","path-length-l","path-length-s (in K-y-space)");

	// Print initialized position for trajectory //
	printf(         "%ld\t%17.17g\t%17.17g\t%g\t%g\n",count_steps,K,y,Deltal,Deltas);
	fprintf(outfile,"%ld\t%17.17g\t%17.17g\t%g\t%g\n",count_steps,K,y,Deltal,Deltas);
	fflush(outfile);

	// Take additional steps //
	if(OpA==1||OpA==3) flowdir = +1;
	if(OpA==2||OpA==4) flowdir = -1;
	dl = flowdir*trajectorydl;
	if(OpA==1||OpA==2){ // const dl
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
				fprintf(outfile,"%ld\t%17.17g\t%17.17g\t%g\t%g\n",count_steps,K,y,Deltal,Deltas);
				printf(          "%ld\t%17.17g\t%17.17g\t%g\t%g\n",count_steps,K,y,Deltal,Deltas);
				fflush(outfile);
			}
		}
	}
	if(OpA==3||OpA==4){ // limited ds
		while(count_steps<Max_steps){ // cutoff arranged in code here !!
			oldK = K;
			oldy = y;
			rk4(EqRecRel, &l, z, dl, N);  // equil: rk4 3D (K,y), EqRecRel
			dlnew = dl;
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
				fprintf(outfile,"%ld\t%17.17g\t%17.17g\t%g\t%g\n",count_steps,K,y,Deltal,Deltas);
				printf(          "%ld\t%17.17g\t%17.17g\t%g\t%g\n",count_steps,K,y,Deltal,Deltas);
				fflush(outfile);
			}
		}
	}

	// Close plotting output file //
	fclose(outfile);

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



// Equilibrium (K,y) recursion relations //
void EqRecRel(double x, double z[2], double dzdx[2], unsigned n){
	dzdx[0] = -A2*z[0]*z[0]*z[1];
	dzdx[1] =  z[1]*(4.0-2.0*PI*z[0]);
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
