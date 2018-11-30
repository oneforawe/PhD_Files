#filename: thesis_plot_fit_Pdep_CcK0c_lmax100_dl0.001_Op2.gnu
reset
set terminal postscript color eps enhanced
set output 'thesis_plot_fit_Pdep_CcK0c_lmax100_dl0.001_Op2.ps'
set size 0.85,0.85

set title "^4He Critical Pair Parameters vs Pressure {/Times-Italic P}"

set ylabel "{/Times-Italic C}_c and {/Times-Italic K}_{0c} (unitless)"
set xlabel "{/Times-Italic P} (bar)"
set xrange [0:30]
#set yrange [0:3e-5]
#set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
set key inside right box width -5 height 0.5 spacing 1.2 title " Vortex loop theory \n calculation fits "
#set rmargin 35

#Old
#Defining the fit parameters and equation:
#C00 = 1.10650545407677
#C01 = -0.0299685857230346
#C02 = 0.000379859757598446
#C03 = -4.09656062944955e-06
#C04 = 1.47249018766356e-07
#C05 = -1.27988117946541e-08
#C06 = 6.43855610603196e-10
#C07 = -1.68170389156538e-11
#C08 = 2.14762244954746e-13
#C09 = -9.71276572314535e-16
#C10 = -1.46448002881887e-18
#fCc(x)  = C00 + C01*x + C02*x**2 + C03*x**3 + C04*x**4 + C05*x**5 + C06*x**6 + C07*x**7 + C08*x**8 + C09*x**9 + C10*x**10

#Defining the fit parameters and equation:
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
fCc(x)  = C00 + C01*x + C02*x**2 + C03*x**3 + C04*x**4 + C05*x**5 + C06*x**6 + C07*x**7 + C08*x**8 + C09*x**9 + C10*x**10

#Defining the fit parameters and equation:
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
fK0c(x)  = k00 + k01*x + k02*x**2 + k03*x**3 + k04*x**4 + k05*x**5 + k06*x**6 + k07*x**7 + k08*x**8 + k09*x**9 + k10*x**10


plot fCc(x)  title " {/Times-Italic C}_c fit ", \
     fK0c(x) title " {/Times-Italic K}_{0c} fit ", \
     'Pdep_K0c_lmax100_dl0.001_Op2.dat' using 4:1 with points linestyle 1 lc rgb "black" notitle, \
     'Pdep_K0c_lmax100_dl0.001_Op2.dat' using 4:2 with points linestyle 1 lc rgb "black" notitle

# Cc points
#     'Pdep_CcK0c.dat' using 1:2 with points linestyle 1 lc rgb "black" notitle, \
# K0c points
#     'Pdep_CcK0c.dat' using 1:3 with points linestyle 1 lc rgb "black" notitle, \
# Cc points
#     'Pdep_K0c_lmax100_dl0.001_Op2.dat' using 4:1 with points linestyle 1 lc rgb "black" notitle, \
# K0c points
#     'Pdep_K0c_lmax100_dl0.001_Op2.dat' using 4:2 with points linestyle 1 lc rgb "black" notitle
# Cc points (redundant)
#     'Pdep_Cc_lmax100_dl0.001_Op2.dat' using 3:1 with points linestyle 1 lc rgb "black" notitle
