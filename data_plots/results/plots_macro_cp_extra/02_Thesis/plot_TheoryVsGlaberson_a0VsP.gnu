#filename: plot_TheoryVsGlaberson_a0VsP.gnu
# This file is copied from where it was in the folder .../programs/results
reset
set terminal postscript color eps enhanced
set output 'plot_TheoryVsGlaberson_a0VsP.ps'

set title "Vortex-Loop Bare Core Diameter {/Times-Italic a}_0({/Times-Italic P}) / {/Times-Italic a}_0(0) (at {/Times-Italic T} = 0.3 K) \n vs Pressure {/Times-Italic P}"

set ylabel "{/Times-Italic a}_0({/Times-Italic P}) / {/Times-Italic a}_0(0) (unitless)"
set xlabel "{/Times-Italic P} (bar)"
set xrange [0:30]
#set yrange [0:0.008]
set grid
set key inside top left box width -4 height 0.5 spacing 1.4
set rmargin 2

PISQ = 9.86960440108935861883

#Parameters & eqn for Tc(P):
# HEY, shouldn't Tc0 be 2.1768 K, for P=0?
Tc0 = 2.17349425585161
Tc1 = -0.00982499579394534
Tc2 = -0.000118194448444384
Tc3 = -4.36914591522034e-07
Tc4 = 7.39407378262721e-09
Tc(y)  = Tc0 + Tc1*y + Tc2*y**2 + Tc3*y**3 + Tc4*y**4

#Parameters & eqn for K0c(P):
k00 = 0.295359928991731
k01 = 0.00508413174676733
k02 = 5.15885451496554e-05
k03 = 6.69186316462475e-07
k04 = -2.74622551950506e-08
k05 = 2.92854219547921e-09
k06 = -8.80012434467914e-11
k07 = -1.93460846518787e-13
k08 = 7.13136752136051e-14
k09 = -1.47121339074112e-15
k10 = 9.60439029245856e-18
K0c(y)  = k00 + k01*y + k02*y**2 + k03*y**3 + k04*y**4 + k05*y**5 + k06*y**6 + k07*y**7 + k08*y**8 + k09*y**9 + k10*y**10

#Parameters & eqn for rho(T,P):
R00 = 145.145109496329
R10 = -0.0976539693059151
R20 = 0.334163407001684
R30 = -0.446930785976304
R40 = 0.181879478545246
R01 = 1.74477604495583
R11 = -0.0919538993179052
R21 = 0.179844560873926
R31 = -0.133606331352667
R41 = 0.0410225514249919
R02 = -0.0491655379690169
R12 = 0.00710698898070406
R22 = -0.00823054225495917
R32 = 0.000609542602247143
R42 = 0.00114916775392305
R03 = 0.0013415037643754
R13 = -0.000362007479155809
R23 = 0.000358809384119286
R33 = 6.48183954357527e-05
R43 = -0.000104112551302631
R04 = -1.69907294147191e-05
R14 = 5.53820368251513e-06
R24 = -3.15773411117433e-06
R34 = -4.99967306908062e-06
R44 = 3.41331223468399e-06
rho(x,y) = R00      + R10*x      + R20*x**2      + R30*x**3      + R40*x**4      + R01*y    + R11*x*y    + R21*x**2*y    + R31*x**3*y    + R41*x**4*y    + R02*y**2 + R12*x*y**2 + R22*x**2*y**2 + R32*x**3*y**2 + R42*x**4*y**2 + R03*y**3 + R13*x*y**3 + R23*x**2*y**3 + R33*x**3*y**3 + R43*x**4*y**3 + R04*y**4 + R14*x*y**4 + R24*x**2*y**4 + R34*x**3*y**4 + R44*x**4*y**4

a0(x,y) = Tc(y)*K0c(y)/rho(x,y)


plot 'Glaberson_a0VsP.dat' using ($1*(300./438)*(1./14.5037738007)):(1+$2*(0.3/350)):(1+$3*(0.3/350)):(1+$4*(0.3/350)) with yerrorbars pt 4 title " Glaberson's data ", \
     '3Dvlt_macro_cp_extra_BareValues_T0.3_Output30.dat' using 1:13 with linespoints title " VLT calculations (T = 0.3 K) ", \
     '3Dvlt_macro_cp_extra_BareValues_Tc_Output30.dat' using 1:13 with linespoints title " VLT calculations (T = Tc) ", \
     a0(0.3,x)/a0(0.3,0) title " Formula using prior fits (T = 0.3 K) ", \
     a0(Tc(x),x)/a0(Tc(x),0) title " Formula using prior fits (T = Tc) "
