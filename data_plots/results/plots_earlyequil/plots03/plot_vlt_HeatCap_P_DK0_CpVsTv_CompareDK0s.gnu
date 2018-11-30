#filename: plot_vlt_HeatCap_P_DK0_CpVsTv_CompareDK0s.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_vlt_HeatCap_P_DK0_CpVsTv_CompareDK0s.ps'

set title "^4He Molar Specific Heat Capacity at Constant Pressure {/Times-Italic c_P} \n vs Temperature Variable {/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l}"

set ylabel "{/Times-Italic c_P} (J mol^{-1} K^{-1})"
set xlabel "{/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l}{/Times-Italic (P)} (unitless)"
set logscale x
set xrange [0.01:1e-6]
#set xrange [1e-2:1e-8]
set yrange [50:200]
set grid
set key inside left top box width -8 spacing 1.2 title "Calculated from \n {/CM-Typewriter vlt\\_HeatCap\\_P.c}"
set rmargin 4

plot 'vlt_HeatCap_P_25.868_DK0_1e-06.dat' using 2:9 with linespoints title " {/Times-Italic P} = 25.868 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.836 K) ", \
     'vlt_HeatCap_P_25.868_DK0_5e-07.dat' using 2:9 with linespoints title " {/Times-Italic P} = 25.868 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.836 K) ", \
     'vlt_HeatCap_P_25.868_DK0_1e-07.dat' using 2:9 with linespoints title " {/Times-Italic P} = 25.868 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.836 K) "
#     -15.8282446947989*log(x) - 16.5178636413462 notitle
