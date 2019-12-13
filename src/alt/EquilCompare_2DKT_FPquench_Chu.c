// EquilCompare_2DKT_FPquench_Chu.c
// This is a version of Chu's program (2DT5_11_00.c), for editing and exploration.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define maxT 5000000	
#define outputT 500
#define outputr 10
#define maxsteps 1000
double const stepsize = 0.008;
double const  PI = 3.141592653589793;
double const  x0=0.0;
double const  k0c = 0.747853;
double const  a0 = 1.0;
double const oneoverTini=1.001001001;
double const qdot = 0.0;				//Turbulence parameter
int i,j,n;
double k[maxsteps+1],g[maxsteps+1];
extern derk();

kt(x,z,zp,n)
double  x,z[2],zp[2];
unsigned n;
{ 	
	zp[0] = -4*PI*PI*PI*z[1]*z[0]*z[0]*exp(4.0*(x));
	// Chu's equil eqn: //
	zp[1] = -qdot-2.0*PI*z[0]*z[1];
	// My new equil eqn: //
//	zp[1] = z[1]*(4.0-4.0*x-2.0*PI*z[0]); // THIS IS WRONG!  I just took the derivative incorrectly to get this!
}

kstepper(x,kstep,kstepp,n)
double x,kstep[1],kstepp[1];
unsigned n;
{ 	
	kstepp[0] = -4*PI*PI*PI*g[j-1]*kstep[0]*kstep[0]*exp(4.0*(x-stepsize));
}

gstepper(x, gin, gstep)
double x, gin[maxsteps+1], gstep[1];
{
	gstep[0] = exp(-2.0*x)*((cDiff/2.0/PI)*(gin[j+1] - 2.0*gin[j] + gin[j-1])/(stepsize*stepsize)
						+(k[j+1]*g[j+1] - k[j-1]*g[j-1])/2.0/stepsize);
}
// NOTE THAT GIN is a modified g, but g is the unmodified g (!)


main (){
	double oneoverT;
	double k0,z[2],x,kstep[1],rhosum,kout;
	double dt,gstep[1];
	double g1[maxsteps+1],g2[maxsteps+1];
	double k1[maxsteps+1],k2[maxsteps+1];
	char rhos[10],gamma[15],r[10],nstring[5];
	FILE *f,*f1;

	// Using Chu's equil eqn (see kt above): //
	f=fopen("EquilCompare_tvdensq2.txt","w");
	f1=fopen("EquilCompare_tgammaq2.txt","w");
	// Using my new equil eqn (see kt above): //
//	f=fopen("EquilCompare_tvdensq2_NewEquilEqn.txt","w");
//	f1=fopen("EquilCompare_tgammaq2_NewEquilEqn.txt","w");

//	Defining the initial condition 	
//
//	alpha = oneoverTini-oneoverTfnl;
	oneoverT = oneoverTini;
	k0 = k0c*oneoverTini;
	z[0]=k0;
	j=0;x=x0;rhosum=0.0;
	k[0] = z[0];
	g[0] = exp(-4.0*x-PI*PI*k0);
	z[1] = g[0];
	rhosum=0.5*g[0]*exp(2.0*x)*stepsize;
	fprintf(f1,"r0    rhos0    gamma0 \n");
	fprintf(f1,"%g %g %g\n",x,k[0]/k0c,g[0]);
	for (j=1;j<maxsteps+1;j++){
	     derk(kt,&x,z,stepsize,2);
	     g[j]=z[1];
	     k[j]=z[0];
	     if(j<maxsteps-4) rhosum=rhosum+exp(2.0*x)*g[j]*stepsize;
	     if(j==maxsteps-4) {rhosum=rhosum+0.5*exp(2.0*x)*g[j]*stepsize;kout=k[j];}
	     if(j%outputr==0) fprintf(f1,"%g %g %g\n",x,k[j]/k0c,g[j]);
	     }
	
	fclose(f);
	fclose(f1);
}
