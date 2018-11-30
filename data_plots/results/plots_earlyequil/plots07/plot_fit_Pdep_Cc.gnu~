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
C00 = 0.255507006991423
C01 = -0.00256258344901229
C02 = -1.53473182251833e-06
C03 = -1.16427278657462e-08
C04 = -4.91527334330677e-11
C05 = -1.82670326116497e-13
C06 = -4.61276030161338e-16
C07 = -6.76444726235525e-19
C08 = -4.26780362659065e-22
fCc(x)  = C00 + C01*x + C02*x**2 + C03*x**3 + C04*x**4 + C05*x**5 + C06*x**6 + C07*x**7 + C08*x**8

plot 'Pdep_Cc.dat' using 3:1 with points, \
     fCc(x)
