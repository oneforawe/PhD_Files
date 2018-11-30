#filename: plot_HeatCapVsTemp.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_HeatCapVsTemp.ps'

set title "^4He Molar Specific Heat Capacity at Constant Volume {/Times-Italic c_P} \n vs Temperature {/Times-Italic T / T}_{/Symbol l}"

set ylabel "{/Times-Italic c_P} (J mol^{-1} K^{-1})"
set xlabel "{/Times-Italic T / T}_{/Symbol l} (unitless)"
set xrange [0.991:1.001]
set yrange [0:150]
set grid
set key inside left top box width -2 height 0.5 title " Calculated using \n {/CM-Typewriter vlt\\_HeatCap.c} "
set rmargin 2

plot 'vlt_HeatCap_DK0_1e-8.dat' using (1-$1):2 title ' DK0 = 1e-8 ' with linespoints, \
     'vlt_HeatCap_DK0_1e-7.dat' using (1-$1):2 title ' DK0 = 1e-7 ' with linespoints, \
     'vlt_HeatCap_DK0_1e-6.dat' using (1-$1):2 title ' DK0 = 1e-6 ' with linespoints, \
     'vlt_HeatCap_DK0_1e-5.dat' using (1-$1):2 title ' DK0 = 1e-5 ' with linespoints, \
     'vlt_HeatCap_DK0_1e-4.dat' using (1-$1):2 title ' DK0 = 1e-4 ' with linespoints, \
     'vlt_HeatCap_DK0_1e-3.dat' using (1-$1):2 title ' DK0 = 1e-3 ' with linespoints, \
     'vlt_HeatCap_DK0_1e-2.dat' using (1-$1):2 title ' DK0 = 1e-2 ' with linespoints, \
     'vlt_HeatCap_DK0_1e-1.dat' using (1-$1):2 title ' DK0 = 1e-1 ' with linespoints
