#filename: plot_vlt_HeatCap_P_00.050_CpVsTvDK0.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_vlt_HeatCap_P_00.050_digits_1_CpVsTvDK0.ps'

set title "^4He Molar Specific Heat Capacity at Constant Pressure {/Times-Italic c_P} \n vs Temperature Variable {/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l} \n using various derivative steps {/Times-Italic DK_0}"

set ylabel "{/Times-Italic c_P} (J mol^{-1} K^{-1})"
set xlabel "{/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l}{/Times-Italic (P)} (unitless)"
set xrange [1:10**(-8)]
set yrange [0:400]
set logscale x
set grid
set key inside left top box height 0.5 spacing 1.1 title " Calculated using \n {/CM-Typewriter vlt\\_HeatCap\\_P.c} \n {/Times-Italic P} = {/Times-Italic P}_{sv} = 0.050 bar "
set rmargin 2

plot 'vlt_HeatCap_P_00.050_DK0_1e-08.dat' using 2:9 title ' DK0 = 1e-8 ' with linespoints, \
     'vlt_HeatCap_P_00.050_DK0_5e-08.dat' using 2:9 title ' DK0 = 5e-8 ' with linespoints, \
     'vlt_HeatCap_P_00.050_DK0_1e-07.dat' using 2:9 title ' DK0 = 1e-7 ' with linespoints, \
     'vlt_HeatCap_P_00.050_DK0_1e-06.dat' using 2:9 title ' DK0 = 1e-6 ' with linespoints, \
     'vlt_HeatCap_P_00.050_DK0_1e-05.dat' using 2:9 title ' DK0 = 1e-5 ' with linespoints, \
     'vlt_HeatCap_P_00.050_DK0_1e-04.dat' using 2:9 title ' DK0 = 1e-4 ' with linespoints, \
     'vlt_HeatCap_P_00.050_DK0_1e-03.dat' using 2:9 title ' DK0 = 1e-3 ' with linespoints

