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

#Defining the fit parameters and equation:
gC00 = 1.08792282377619
gC01 = -0.0292860526555516
gC02 = 0.000366954878475542
gC03 = -3.77409045188115e-06
gC04 = 1.36370669649129e-07
gC05 = -1.38573433201423e-08
gC06 = 8.00276663858592e-10
gC07 = -2.47619512902066e-11
gC08 = 4.15339417790346e-13
gC09 = -3.49992602072889e-15
gC10 = 1.1221557142682e-17
gCc(x)  = gC00 + gC01*x + gC02*x**2 + gC03*x**3 + gC04*x**4 + gC05*x**5 + gC06*x**6 + gC07*x**7 + gC08*x**8 + gC09*x**9 + gC10*x**10


plot 'Pdep_Cc_lmax100_dl0.001_Op2.dat' using 3:1 with points title " Calculated points ", \
     gCc(x) title " New Fit ", \
     fCc(x) title " Old Fit "
