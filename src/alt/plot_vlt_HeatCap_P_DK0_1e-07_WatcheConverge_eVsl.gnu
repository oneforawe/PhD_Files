#filename: plot_vlt_HeatCap_P_DK0_1e-07_WatcheConverge_eVsl.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_vlt_HeatCap_P_DK0_1e-07_WatcheConverge_eVsl.ps'

set title "^4He Helmholtz Parameter {/Times-Italic e} vs Length Scale {/Times-Italic l} \n Watching {/Times-Italic e} Converge in Runge Kutta Loop"

set ylabel "{/Times-Italic e} (unitless)"
set xlabel "{/Times-Italic l} (unitless)"
#set logscale x
#set xrange [0.75:0.95]
#set yrange [-0.08:0]
set grid
set key inside right top box width -4 spacing 1.2 title "Calculated from \n {/CM-Typewriter vlt\\_HeatCap\\_P.c}"
set rmargin 2

plot 'vlt_HeatCap_P_00.050_DK0_1e-07_WatcheConverge.dat' using 3:6 with linespoints title " {/Times-Italic P} = {/Times-Italic P}_{sv} = 0.050 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.172 K) "
