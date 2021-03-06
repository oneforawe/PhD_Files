#filename: plot_fit_Pdep_Cc.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_fit_Pdep_Cc.ps'
#set size 1.0,1.1

set title "^4He Core parameter {/Times-Italic C} vs Pressure {/Times-Italic P}"

set ylabel "{/Times-Italic C} (unitless)"
set xlabel "{/Times-Italic P} (bar)"
#set xrange [0:4.38e-4]
#set yrange [0:3e-5]
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
set key inside right box width 0 height 0.5 spacing 1.2 title " VLT calculations \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.67168835 "
#set rmargin 35

#Defining the fit parameters and equation:
C00 = 1.10650545407677
C01 = -0.0299685857230346
C02 = 0.000379859757598446
C03 = -4.09656062944955e-06
C04 = 1.47249018766356e-07
C05 = -1.27988117946541e-08
C06 = 6.43855610603196e-10
C07 = -1.68170389156538e-11
C08 = 2.14762244954746e-13
C09 = -9.71276572314535e-16
C10 = -1.46448002881887e-18
fCc(x)  = C00 + C01*x + C02*x**2 + C03*x**3 + C04*x**4 + C05*x**5 + C06*x**6 + C07*x**7 + C08*x**8 + C09*x**9 + C10*x**10

plot 'Pdep_Cc.dat' using 3:1 with points title " Calculated points ", \
     fCc(x) title " Fit "
