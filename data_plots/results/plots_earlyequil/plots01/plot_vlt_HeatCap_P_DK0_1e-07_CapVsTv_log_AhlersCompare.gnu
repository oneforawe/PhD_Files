#filename: plot_vlt_HeatCap_P_DK0_1e-07_CapVsTv_log_AhlersCompare.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_vlt_HeatCap_P_DK0_1e-07_CapVsTv_log_AhlersCompare.ps'

set title "^4He Molar Specific Heat Capacity (at const P or V?) {/Times-Italic c} \n vs Temperature Variable {/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l}"

set ylabel "{/Times-Italic c} (units not right yet)"
set xlabel "{/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l}({/Times-Italic P}) (unitless)"
set logscale x
set xrange [1:1e-7]
#set xrange [1e-2:1e-8]
set yrange [0:120]
set grid
set key inside left top box width -4 spacing 1.2 title "Calculated from \n {/CM-Typewriter vlt\\_HeatCap\\_P.c}"
set rmargin 2

plot 'vlt_HeatCap_P_00.050bar_DK0_1e-07.dat' using 1:2 with linespoints title " {/Times-Italic P} = {/Times-Italic P}_{sv} = 0.050 bar ", \
     'vlt_HeatCap_P_01.646bar_DK0_1e-07.dat' using 1:2 with linespoints title " {/Times-Italic P} = 1.644 bar ", \
     'vlt_HeatCap_P_07.328bar_DK0_1e-07.dat' using 1:2 with linespoints title " {/Times-Italic P} = 7.328 bar ", \
     'vlt_HeatCap_P_15.031bar_DK0_1e-07.dat' using 1:2 with linespoints title " {/Times-Italic P} = 15.031 bar ", \
     'vlt_HeatCap_P_18.180bar_DK0_1e-07.dat' using 1:2 with linespoints title " {/Times-Italic P} = 18.180 bar ", \
     'vlt_HeatCap_P_22.533bar_DK0_1e-07.dat' using 1:2 with linespoints title " {/Times-Italic P} = 22.533 bar ", \
     'vlt_HeatCap_P_25.868bar_DK0_1e-07.dat' using 1:2 with linespoints title " {/Times-Italic P} = 25.868 bar "
