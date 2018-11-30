#filename: plot_fit_Pdep_Cc.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_fit_Pdep_Cc.ps'
set size 1.0,1.1

set title "^4He Core quantity {/Times-Italic C}_c vs Pressure {/Times-Italic P}"

set ylabel "{/Times-Italic C}_c (unitless)"
set xlabel "{/Times-Italic P} (bar)"
#set xrange [0:4.38e-4]
#set yrange [0:3e-5]
set grid
set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
set rmargin 35

#Defining the fit parameters and equation:
C00 = 1.10173855010172
C01 = -0.0280346715561482
C02 = 0.000335204706379962
C03 = -3.22365242538447e-06
C04 = -1.94786927555144e-08
C05 = 1.01992137142591e-08
C06 = -1.11673138328809e-09
C07 = 7.15793487126364e-11
C08 = -2.90383266833513e-12
C09 = 7.53075469261995e-14
C10 = -1.20907488333841e-15
C11 = 1.09300580803761e-17
C12 = -4.24495149368875e-20
fCc(x)  = C00 + C01*x + C02*x**2 + C03*x**3 + C04*x**4 + C05*x**5 + C06*x**6 + C07*x**7 + C08*x**8 + C09*x**9 + C10*x**10 + C11*x**11 + C12*x**12

plot 'Pdep_Cc.dat' using 3:1 with points, \
     fCc(x)
