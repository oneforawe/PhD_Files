/*  File Comments  */
/*=========================================================================*/

#include <stdio.h>
#define IN    // input label
#define OUT   // output label
#define INOUT // input/output label
#define N     2  // If using ring1.c N=3, If using 2DKT.c or K0cFind1.c N=2

/* Data type definitions */
typedef struct {
	double c[N];
} RETARRAY;



/*  Function definition  */
/*=========================================================================*/
//   derk1(func(),calculated,Dl,n)
//   func(l,z)

double derk1(double (*func)(double l, double z[]), RETARRAY calcKG, double calcl, double Dl, unsigned n){
	int i;
	double dl, testl[3], DKDl[4];
	RETARRAY testK[3];

	/* Calculate the derivatives at the input calculated state and 3n test states */
	DKDl[0] = func(IN calcl, IN calcKG.c);
	for(i=0; i<3; i++){
		dl = Dl*((i+2)/2)/2;  // See "Note" below.
		testK[i].c[0] = calcKG.c[0] + DKDl[i]*dl;
		testK[i].c[1] = calcKG.c[1];
		testl[i] = calcl + Dl*((i+1)/2)/2;  // See "Note" below.
		DKDl[i+1] = func(IN testl[i], IN testK[i].c);
	}
	// Note: This algorithm, using '((...)/2)/2', creates a desired rounding "error" that produces the sequence 0.0, 0.5, 0.5, 1.0, instead of 0.0, 0.5, 0.75, 1.0.

	/* Update the calculated state, using the derivatives calculated above */
	calcKG.c[0] += Dl*(DKDl[0]+2*(DKDl[1]+DKDl[2])+DKDl[3])/6;
	// calculated.l += Dl;

	/* Return the new calculated K-value (not whole state) */
	return calcKG.c[0];
}



/*  Program Notes  */
/*=========================================================================*/
/*

See further documentation (UCLAresearchNotes.pdf, UCLAresearchPrograms.pdf, vlt_ThermStates.c, ring1.c) for elaboration.

Dl is the increment in l.
n is the number of thermodynamic variables (n=3 given thermv[i], i=1,2,3).

DKDl[i].arr[0] = dK/dl
DKDl[i].arr[1] = dy/dl
DKDl[i].arr[2] = de/dl

*/
