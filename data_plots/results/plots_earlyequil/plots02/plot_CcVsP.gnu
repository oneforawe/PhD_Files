#filename: plot_CcVsP.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_CcVsP.ps'

set title "^4He Critical Vortex-Core Parameter {/Times-Italic C}_c \n vs Pressure {/Times-Italic P}"

set ylabel "{/Times-Italic C}_c (unitless)"
set xlabel "{/Times-Italic P} (bar)"
#set xrange [1:1e-8]
#set mxtics 10
#set yrange [0:300]
set grid
set key inside right top box width -6 spacing 1.2 title "Curve and points from \n {/CM-Typewriter vlt\\_K0cFind.c} and \n Greywall/Ahlers fits"
set rmargin 4

#From fit_PressureDependences_2.par
Cc00 =   1.10652102266763
Cc01 =  -0.0299766068265842
Cc02 =   0.000380954449384652
Cc03 =  -3.60465470172711e-06
Cc04 =  -1.13928690475838e-08
Cc05 =   5.27191245296858e-09
Cc06 =  -4.02970077711438e-10
Cc07 =   1.8784597083754e-11
Cc08 =  -5.69181439968371e-13
Cc09 =   1.12069994170299e-14
Cc10 =  -1.38157107961357e-16
Cc11 =   9.67488197927817e-19
Cc12 =  -2.93300475897348e-21
fCc(x) = Cc00 + Cc01*x + Cc02*x**2 + Cc03*x**3 + Cc04*x**4 + Cc05*x**5 + Cc06*x**6 + Cc07*x**7 + Cc08*x**8 + Cc09*x**9 + Cc10*x**10 + Cc11*x**11 + Cc12*x**12

plot fCc(x) title " 12th-order fit ", \
     'Pdep.dat' using 1:3 title " Points calc'd from fit "