// 2DKT_FPquench_simplec.c

#include <stdio.h>
#include <math.h>

const double PI   = 3.14159265358979323846;
const double PISQ = 9.86960440108935799230;
const double PICU = 31.00627668029981620634;
const double B    = 4.0*PICU;
const int    N    = 2;

const int    lmax     = 10;          //  1
const int    lsteps   = 100000;      //  lsteps = lmax/dl = 10/0.0001 = 100,000  (make sure this is divisible by ldatamax)
const int    lpts     = lsteps+1;    //  from l=0 to l=lmax, inclusive (incl. the boundary-condition-enforcing pnt, there are lpts+1 "l" pnts)
const int    ldatamax = 100;         //  max number of data points recorded per time step
// const double dt    = 1;           //  time step (in units of "diffusion time"(?))
const int    tmax     = 5;           //  (later, 10000) max unitless time (see if(t==...) below)
// tdatamax           = 500;?
const double iTfrac   = 0.95;        //  initial temperature fraction Ti/Tkt
const double qTfrac   = 0.10;        //  quench temperature fraction Tq/Tkt
const double a0       = 1.0;
const double a04      = 1.0;         //  a0 to the fourth power
const double K0c      = 0.747853;

int i,j,n;
double K[lsteps+1],G[lsteps+1];


/*  Function prototypes  */
/*=========================================================================*/
//  derk(f,l,y,h,n)
//  Kstepper(l,Kstep,Kstepp,n)
//  Gstepper(l,Gin,Gstep)
//  Kt(l,z,zp,n)
int derk(void (*f)(), double *l, double y[], double h, unsigned n);
void Kstepper(double l, double Kstep[1], double Kstepp[1], unsigned n);
void Gstepper(double l, double Gin[lsteps+1], double Gstep[1]);
void Kt(double l, double z[2], double zp[2], unsigned n);


main(){
	double K0,y,z[2],l,Kstep[1],rhosum,Kout;
	double t, dt, dtLimit;
	double const l0=0.0;
	double Gstep[1];
	double G1[lsteps+1],G2[lsteps+1],G3[lsteps+1];
	double k1[lsteps+1],k2[lsteps+1],k3[lsteps+1],k4[lsteps+1];

	FILE *outfile_ProbDensG, *outfile_VPairDens;
	outfile_ProbDensG = fopen("2DKT_FPquench_simpleb_ProbDensG.dat","w");
	outfile_VPairDens = fopen("2DKT_FPquench_simpleb_VPairDens.dat","w");
	fprintf(outfile_ProbDensG, "l    rhos    Gamma\n");
	printf(                    "l    rhos    Gamma\n");
	fprintf(outfile_VPairDens, "?   rhosum  Kout/K0c\n");


	// Defining the initial condition

	t=0.0;
	l=l0;rhosum=0.0;
	K0 = K0c/iTfrac;
	y=exp(-PISQ*K0/2);
	z[0]=K0;
	K[0]=z[0];
	z[1]=y;
	G[0]=z[1]*z[1]*exp(-4.0*l);

	fprintf(f1,"%g\t%g\t%g\n",l,K[0]/K0,G[0]);

	for(j=1; j<lsteps+1; j++){
		derk(Kt,&l,z,dl,2);
		G[j]=z[1]*z[1]*exp(-4.0*l);
		K[j]=z[0]; 
		if(j%(lsteps/datamax)==0||j<10) // (j%100==0)
			fprintf(f1,"%g\t%g\t%g\n",l,K[j]/K0,G[j]);
			printf(    "%g\t%g\t%g\n",l,K[j]/K0,G[j]);
	}
	fprintf(f,"%g\t%g\t%g\n",0*Kout,rhosum,Kout/K0c);
	printf(   "%g\t%g\t%g\n",0*Kout,rhosum,Kout/K0c);


	// Quench to qTfrac (not oneoverT)

	l=l0;
	K[0]=K0c/qTfrac;
	Kstep[0]=K[0];
	for(j=1; j<lsteps+1; j++){
		derk(Kstepper,&l,Kstep,dl,1);
		K[j]=Kstep[0];
	}

	fprintf(f1,"\n%s\n","First adjustment of K/K0:");
	l=l0;
	for(j=0; j<lsteps+1; j++){
		if(j%(lsteps/datamax)==0||j<10)
			fprintf(f1,"%g\t%g\t%g\n",l,K[j]/K[0],G[j]);
		l += dl;
	}


	for(j=0;j<=lsteps;j++){
		G1[j] = 0.0;
		G2[j] = 0.0;
		G3[j] = 0.0;
	}

	/*  End of defining the initial condition  */


	// Define the timestep size

	dtLimit = 2.5*(dl*dl)/sqrt((4.0/PI/PI) + K[0]*dl*K[0]*dl);
	dt = 0.9*dtLimit;


	// Step out in time

	for(n=1;n<=tmax;n++){
		t += dt;
		fprintf(f1,"\n# n = %i\tt = %f\n\n", n, t);  /* */

		l=l0+dl;
		for(j=1; j<lsteps; j++){
			Gstepper(l,G,Gstep);
			k1[j] = Gstep[0]*dt;
			G1[j] = G[j] + k1[j]/2.0;
			l += dl;
		}
		G1[0] = G1[1];
		G1[lsteps] = G1[lsteps-1];

		l=l0+dl;
		for(j=1; j<lsteps; j++){
			Gstepper(l,G1,Gstep);
			k2[j] = Gstep[0]*dt;
			G2[j] = G[j] - k1[j] + k2[j];
			l += dl;
		}
		G2[0] = G2[1];
		G2[lsteps] = G2[lsteps-1];

		l=l0+dl;
		for(j=1; j<lsteps; j++){
			Gstepper(l,G2,Gstep);
			k3[j] = Gstep[0]*dt;
//			G3[j] = G[j] + k1[j]/6.0 + 2.0*k2[j]/3.0 + k3[j]/6.0;
			G[j] = G[j] + k1[j]/6.0 + 2.0*k2[j]/3.0 + k3[j]/6.0;
			l += dl;
		}
/*		G3[0] = G3[1];
		G3[lsteps]= G3[lsteps-1];

		l=l0;
		for(j=1;j<lsteps;j++){
			Gstepper(l,G3,Gstep);
			k4[j] = Gstep[0]*dt;
			G[j] = G[j] + k1[j]/6.0 + k2[j]/3.0 + k3[j]/3.0 + k4[j]/6.0;
			l += dl;
		}  */
		G[0] = G[1];
		G[lsteps] = G[lsteps-1];


		fprintf(f1,"\n%s\n","After adjusting G:");
		l=l0;
		for(j=0; j<lsteps+1; j++){
			if(j%(lsteps/datamax)==0||j<10)
				fprintf(f1,"%g\t%g\t%g\n",l,K[j]/K[0],G[j]);
			l += dl;
		}


		// Find K by derk and
		// Output the results

		l=l0;
		K[0]=K0c/qTfrac;
		Kstep[0]=K[0];
		fprintf(f1,"%g\t%g\t%g\n",l,K[0]/K[0],G[0]);

		for(i=1;i<=lsteps;i++){
			derk(Kstepper,&l,Kstep,dl,1);
			K[i]=Kstep[0];
			if(i%(lsteps/datamax)==0||i<10)
				fprintf(f1,"%g\t%g\t%g\n",l,K[i]/K[0],G[i]);
		}

		fflush(f);
		fflush(f1);

	}

	fclose(f);
	fclose(f1);
}



int derk(void (*f)(), double *l, double y[], double h, unsigned n){
	int i,j;
	double *k[4],*s,temp;

	for(i=0; i<4; i++)
	  k[i] = (double *)calloc(n,sizeof(double));
	s = (double *)calloc(n,sizeof(double));
	f(*l,y,k[0],n);
	for(i=1; i<4; i++){
		temp = h*((i+1)/2)/2;
		for(j=0; j<n; j++)
		  s[j] = y[j]+k[i-1][j]*temp;
		f((*l)+h*(i/2)/2,s,k[i],n);
	}
	for(j=0; j<n; j++)
	  y[j] += h*(k[0][j]+2*(k[1][j]+k[2][j])+k[3][j])/6;
	*l += h;
	for(i=0; i<4; i++)
	  free((char *) k[i]);
	free((char *)s);
	return(0);
}



void Kstepper(double l, double Kstep[1], double Kstepp[1], unsigned n){
	Kstepp[0] = -4*PI*PI*PI*G[j-1]*Kstep[0]*Kstep[0]*exp(4.0*l);
}



void Gstepper(double l, double Gin[lsteps+1], double Gstep[1]){
	Gstep[0] = exp(-2.0*l)*((1/2.0/PI)*(Gin[j+1] - 2.0*Gin[j] + Gin[j-1])/(dl*dl) + (K[j+1]*G[j+1] - K[j-1]*G[j-1])/2.0/dl);
}



void EqRecRel(double l, double z[2], double zp[2], unsigned n){
	zp[0] = -4*PI*PI*PI*z[1]*z[1]*z[0]*z[0];
	zp[1] = (2-PI*z[0])*z[1];
}
