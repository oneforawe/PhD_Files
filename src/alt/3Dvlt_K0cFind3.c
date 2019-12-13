//  File Comments  //
//=========================================================================//

/* FILENAME: 3Dvlt_K0cFind3.c */
/* VERSION: 1 (2012 Jun 07 - ...)
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
     and then type "g++ -lm vlt_K0cFind3.c" without the quotes; then to run, type "./a.out".
*/
/* PROGRAM IDEAS:
   * To make this code more efficient, realize that there is really only one trajectory that determines all of the K0c.  The list of Cc's can be loaded all at once, arranged into two lists (one for the positive-K-direction trajectory portion, and the other for the negative-K-direction trajectory portion), and each K0c can be calculated along the way as the trajectory is traversed.
*/

/* RESULTS:
*/



//  Function Preparation  //
