#filename: splot_fit_BrooksDonnelly_RhoVsTP_4order.gnu
reset
set terminal postscript color eps enhanced
set output 'splot_fit_BrooksDonnelly_RhoVsTP_4order.ps'
#set size 1.3,1

set title "^4He Density {/Symbol-Oblique r} \n vs Temperature {/Times-Italic T} and Pressure {/Times-Italic P} "

set zlabel "{/Symbol-Oblique r} (kg m^{-3})"
set xlabel "{/Times-Italic P} (bar)"
set ylabel "{/Times-Italic T} (K)"
#set logscale x
set xrange [0:2.5]
set yrange [0:30]
set zrange [130:180]
#set grid
#set key top box spacing 1.2 title " Data from Brooks/Donnelly "
set key top box spacing 1.2 title " Fit of data from Brooks/Donnelly "
set lmargin 8

set pm3d
splot 'fit_BrooksDonnelly_RhoVsTP_4order_3d.dat' with pm3d


#R00 = 145.145109496329000
#R10 = -0.097653969305915
#R20 = 0.334163407001684
#R30 = -0.446930785976304
#R40 = 0.181879478545246
#R01 = 1.744776044955830
#R11 = -0.091953899317905
#R21 = 0.179844560873926
#R31 = -0.133606331352667
#R41 = 0.041022551424992
#R02 = -0.049165537969017
#R12 = 0.007106988980704
#R22 = -0.008230542254959
#R32 = 0.000609542602247
#R42 = 0.001149167753923
#R03 = 0.001341503764375
#R13 = -0.000362007479156
#R23 = 0.000358809384119
#R33 = 0.000064818395436
#R43 = -0.000104112551303
#R04 = -0.000016990729415
#R14 = 0.000005538203683
#R24 = -0.000003157734111
#R34 = -0.000004999673069
#R44 = 0.000003413312235


#rho(x,y) =  #R00      + #R10*x      + #R20*x**2      + #R30*x**3      + #R40*x**4      + #R01*y    + #R11*x*y    + #R21*x**2*y    + #R31*x**3*y    + #R41*x**4*y    + #R02*y**2 + #R12*x*y**2 + #R22*x**2*y**2 + #R32*x**3*y**2 + #R42*x**4*y**2 + #R03*y**3 + #R13*x*y**3 + #R23*x**2*y**3 + #R33*x**3*y**3 + #R43*x**4*y**3 + #R04*y**4 + #R14*x*y**4 + #R24*x**2*y**4 + #R34*x**3*y**4 + #R44*x**4*y**4 

#set pm3d
#splot rho(x,y) with pm3d
