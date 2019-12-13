#include <stdio.h>
#include <math.h>

#define PI        3.14159265358979323846
#define PISQ      9.86960440108935861883
#define A         4.0*PI*PI*PI/3.0
#define Cc        1.03
#define precision 1e-14
#define K01       0.3

typedef struct {
	double l;
	double tempv;
	double thrmv[2];
} STATE;
typedef struct {
	double arr[2];
} RETARRAY;

extern STATE derk(RETARRAY (*func)(double l, double z[]), STATE calculated, double Dl, unsigned n);
RETARRAY ringrecrel(double l, double z[]);

main(){
	double K0,DK0,oldK,oldy;
	STATE estimated;

	K0 = K01;
	DK0 = 0.1;

	while(DK0>precision){

		estimated.l = 0.0;

		estimated.thrmv[0] = K0;
		estimated.thrmv[1] = exp(-PISQ*K0*Cc);
		oldK = estimated.thrmv[0];
		oldy = estimated.thrmv[1];

		while(1==1){
			estimated = derk(ringrecrel, estimated, 0.0001, 2);

			if(estimated.thrmv[0]>5.0) {
				printf("K exploded\n");
				DK0 = DK0/2;
				K0 = K0-DK0;
				break;
			}
			if(estimated.l>6.0 && estimated.thrmv[0]<oldK) {
				printf("K decreasing\n");
				DK0 = DK0/2;
				K0 = K0+DK0;
				break;
			}
			/* Both of the following if statements should be left out. */
			/* equivalent to Gary's program (redundant with K decreasing?)
			if(estimated.l>6.0 && estimated.thrmv[1]<oldy) {
				printf("y decreasing\n");
				DK0 = DK0/2;
				K0 = K0-DK0;
				break;
			} */
			/* what I would expect to do... (but it's wrong since K may increase before going to zero)
			if(estimated.l>6.0 && estimated.thrmv[0]>oldK) {
				printf("K increasing\n");
				DK0 = DK0/2;
				K0 = K0-DK0;
				break;
			} */
			oldK = estimated.thrmv[0];
			oldy = estimated.thrmv[1];
		}
		printf("%1.15e %1.15e\n",K0,DK0);
	}
}

RETARRAY ringrecrel(double l, double z[]){
	RETARRAY dzdl;
	dzdl.arr[0] = z[0]-A*z[1]*z[0]*z[0];
	dzdl.arr[1] = z[1]*(6.0-PISQ*z[0]*(1.0-0.6*log(z[0])));
	return dzdl;
}
