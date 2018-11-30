#filename: plot_vlt_HeatCap_P_DK0_1e-07_CapVsTv.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_vlt_HeatCap_P_DK0_1e-07_CapVsTv.ps'

set title "^4He Molar Heat Capacity (at const P or V?) {/Times-Italic c} \n vs Temperature Variable {/Symbol-Oblique t}"

set ylabel "{/Times-Italic c} (units not right yet)"
set xlabel "{/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l}({/Times-Italic P}) (unitless)"
#set logscale x
#set xrange [1e-2:1e-8]
#set yrange [0:0.01]
set grid
set key inside right box width -4 spacing 1.2 title " Calculated from \n {/CM-Typewriter vlt\\_HeatCap\\_P.c}"
set rmargin 2

plot 'vlt_HeatCap_P_00.00bar_DK0_1e-07.dat' using 1:2 with linespoints title " {/Times-Italic P} = 00.00 bar ", \
     'vlt_HeatCap_P_00.05bar_DK0_1e-07.dat' using 1:2 with linespoints title " {/Times-Italic P} = {/Times-Italic P}_{sv} = 00.05 bar ", \
     'vlt_HeatCap_P_00.30bar_DK0_1e-07.dat' using 1:2 with linespoints title " {/Times-Italic P} = 00.30 bar ", \
     'vlt_HeatCap_P_00.60bar_DK0_1e-07.dat' using 1:2 with linespoints title " {/Times-Italic P} = 00.60 bar ", \
     'vlt_HeatCap_P_00.70bar_DK0_1e-07.dat' using 1:2 with linespoints title " {/Times-Italic P} = 00.70 bar ", \
     'vlt_HeatCap_P_00.71bar_DK0_1e-07.dat' using 1:2 with linespoints title " {/Times-Italic P} = 00.71 bar ", \
     'vlt_HeatCap_P_08.90bar_DK0_1e-07.dat' using 1:2 with linespoints title " {/Times-Italic P} = 08.90 bar ", \
     'vlt_HeatCap_P_09.00bar_DK0_1e-07.dat' using 1:2 with linespoints title " {/Times-Italic P} = 09.00 bar ", \
     'vlt_HeatCap_P_10.00bar_DK0_1e-07.dat' using 1:2 with linespoints title " {/Times-Italic P} = 10.70 bar ", \
     'vlt_HeatCap_P_15.00bar_DK0_1e-07.dat' using 1:2 with linespoints title " {/Times-Italic P} = 15.00 bar ", \
     'vlt_HeatCap_P_20.00bar_DK0_1e-07.dat' using 1:2 with linespoints title " {/Times-Italic P} = 20.00 bar ", \
     'vlt_HeatCap_P_24.80bar_DK0_1e-07.dat' using 1:2 with linespoints title " {/Times-Italic P} = 24.80 bar "
