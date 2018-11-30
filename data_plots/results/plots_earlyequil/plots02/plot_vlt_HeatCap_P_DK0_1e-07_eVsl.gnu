#filename: plot_vlt_HeatCap_P_DK0_1e-07_eVsl.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_vlt_HeatCap_P_DK0_1e-07_eVsl.ps'

set title "^4He Helmholtz Parameter {/Times-Italic e} vs Length Scale {/Times-Italic l}"

set ylabel "{/Times-Italic e} (unitless)"
set xlabel "{/Times-Italic l} (unitless)"
#set logscale x
#set xrange [0.75:0.95]
set yrange [-0.08:0]
set grid
set key inside right top box width -4 spacing 1.2 title "Calculated from \n {/CM-Typewriter vlt\\_HeatCap\\_P.c}"
set rmargin 2

plot 'vlt_HeatCap_P_00.050_DK0_1e-07.dat' using 3:5 with linespoints title " {/Times-Italic P} = {/Times-Italic P}_{sv} = 0.050 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.172 K) ", \
     'vlt_HeatCap_P_01.646_DK0_1e-07.dat' using 3:5 with linespoints title " {/Times-Italic P} = 1.644 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.157 K) ", \
     'vlt_HeatCap_P_07.328_DK0_1e-07.dat' using 3:5 with linespoints title " {/Times-Italic P} = 7.328 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.095 K) ", \
     'vlt_HeatCap_P_15.031_DK0_1e-07.dat' using 3:5 with linespoints title " {/Times-Italic P} = 15.031 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.998 K) ", \
     'vlt_HeatCap_P_18.180_DK0_1e-07.dat' using 3:5 with linespoints title " {/Times-Italic P} = 18.180 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.954 K) ", \
     'vlt_HeatCap_P_22.533_DK0_1e-07.dat' using 3:5 with linespoints title " {/Times-Italic P} = 22.533 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.889 K) ", \
     'vlt_HeatCap_P_25.868_DK0_1e-07.dat' using 3:5 with linespoints title " {/Times-Italic P} = 25.868 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.836 K) "
