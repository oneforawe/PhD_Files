#include <stdio.h>#include <math.h>#define PI	3.14159265#define PISQ  9.8696044#define A	4.0*PI*PI*PI/3.0#define k0c  0.3091469984216758#define C 1.03#define delt 0.00000003#define a0 2.56extern derk();double k0,a,ac;ring(t,y,yp,n)double t,y[3],yp[3];unsigned n;{	yp[0]=y[0]-y[0]*y[0]*A*y[1];	yp[1]=y[1]*(6.0-PISQ*y[0]*(1.0-0.6*log(y[0])));	yp[2]= -PI*y[1]*exp(-3.0*t);}main()	{	double t,y[3],yp[3],kr,oldkr,yplus[3],ypplus[3],yminus[3],ypminus[3],	k0plus,k0minus,tplus,tminus,cap,temp,D,a,derivplus,derivminus;	FILE *f;	int i,j,k;	f=fopen("caphe.out","a");	for(i=0;i < 35;i++){		temp=1.0e-8*exp(0.4*i);		k0=k0c/(1.0-temp);		a=a0;		y[0]=k0;		yplus[0]=k0+delt;		yminus[0]=k0-delt;		y[1]=exp(-k0*PISQ*C);		yplus[1]=exp(-(k0+delt)*PISQ*C);		yminus[1]=exp(-(k0-delt)*PISQ*C);		y[2]=0.0;		yplus[2]=0.0;		yminus[2]=0.0;		D=-1.0;		j=0;		t=0.0;		tplus=0.0;		tminus=0.0;		kr=k0;		oldkr=k0+1.0;		while ( a>0.0) {			derk(ring,&t,y,0.0001,3);			derk(ring,&tplus,yplus,0.0001,3);			derk(ring,&tminus,yminus,0.0001,3);			a=a0*exp(t);			kr=y[0]*a0/a;			j = j+1;			derivplus=(yplus[2]-y[2])/delt;			derivminus=(y[2]-yminus[2])/delt;			cap=D*k0*k0*(derivplus-derivminus)/delt;			/*if(j%2000==0)				printf("%G %g %g %g \n%1.18e\n",t,kr/k0,y[0],y[1],cap);			if(j>10000)				j=0;*/			if(y[0]>6)				break;		}		cap=D*k0*k0*(derivplus-derivminus)/delt;		fprintf(f,"%G %G %G %G %G\n",temp,cap*8.31,t,kr/k0,y[2]);		printf("%G %G %G %G %G\n",temp,cap*8.31,t,kr/k0,y[2]);		fflush(f);	}	fclose(f);}