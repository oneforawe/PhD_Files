#filename: thesis_plot_TheoryVsBrooksDonnelly_U0VsP.gnu
#This file/plot is similar to plot_RotonEnergyAtTcVsPressure.png and plot_RotonEnergyVsTemp_atTc.png in the folder .../programs/results
reset
set terminal postscript color eps enhanced
set output 'thesis_plot_TheoryVsBrooksDonnelly_U0VsP.ps'
set size 0.85,0.85

set title "^4He Crititcal (near {/Times-Italic T}_{/Symbol-Oblique l}) Smallest Vortex-Loop Energy {/Times-Italic U}_0 (in kelvins) vs Pressure {/Times-Italic P}"

set ylabel "{/Times-Italic U}_0/{/Times-Italic k}_B (K)"
set xlabel "{/Times-Italic P} (bar)"
set xrange [0:30]
#set yrange [0:0.008]
#set grid
set key inside top right box width -12 height 0.5 spacing 1.4
set rmargin 2

PISQ = 9.86960440108935861883

#Parameters & eqn for Tc(P):
Tc0 = 2.17349425585161
Tc1 = -0.00982499579394534
Tc2 = -0.000118194448444384
Tc3 = -4.36914591522034e-07
Tc4 = 7.39407378262721e-09
Tc(x)  = Tc0 + Tc1*x + Tc2*x**2 + Tc3*x**3 + Tc4*x**4

#Parameters & eqn for Cc(P):
C00 = 1.10632942839817
C01 = -0.0299615790448538
C02 = 0.000379514388508258
C03 = -4.10119230465265e-06
C04 = 1.58501290197793e-07
C05 = -1.45151069153968e-08
C06 = 7.65278104409589e-10
C07 = -2.15547139971012e-11
C08 = 3.19404125737557e-13
C09 = -2.19567729108798e-15
C10 = 4.42335546044107e-18
Cc(x)  = C00 + C01*x + C02*x**2 + C03*x**3 + C04*x**4 + C05*x**5 + C06*x**6 + C07*x**7 + C08*x**8 + C09*x**9 + C10*x**10

#Parameters & eqn for K0c(P):
k00 = 0.295390109125547
k01 = 0.00508440693337029
k02 = 5.15936104537981e-05
k03 = 6.6840260137999e-07
k04 = -2.73682424131871e-08
k05 = 2.92445414526824e-09
k06 = -8.79753799932285e-11
k07 = -1.89938965286979e-13
k08 = 7.1196156273502e-14
k09 = -1.46975149440742e-15
k10 = 9.59795096881752e-18
K0c(x)  = k00 + k01*x + k02*x**2 + k03*x**3 + k04*x**4 + k05*x**5 + k06*x**6 + k07*x**7 + k08*x**8 + k09*x**9 + k10*x**10

U0ok(x) = PISQ*K0c(x)*Cc(x)*Tc(x)

plot 'BrooksDonnelly_RotonEnergyVsPressTc.dat' using 1:3 title " Roton data taken/extrapolated from Brooks, Donnelly" with linespoints, \
     U0ok(x) title " Loop theory fit ", \
     '3Dvlt_macro_cp_extra_BareValues_T0.3_OutputAll.dat' using 1:11 with points linestyle 1 lc rgb "black" notitle
