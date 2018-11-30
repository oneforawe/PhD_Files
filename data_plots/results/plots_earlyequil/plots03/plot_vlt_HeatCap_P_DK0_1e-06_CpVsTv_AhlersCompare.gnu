#filename: plot_vlt_HeatCap_P_DK0_1e-06_CpVsTv_AhlersCompare.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_vlt_HeatCap_P_DK0_1e-06_CpVsTv_AhlersCompare.ps'

set title "^4He Molar Specific Heat Capacity at Constant Pressure {/Times-Italic c_P} \n vs Temperature Variable {/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l}"

set ylabel "{/Times-Italic c_P} (J mol^{-1} K^{-1})"
set xlabel "{/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l}{/Times-Italic (P)} (unitless)"
set logscale x
set xrange [1:1e-8]
#set xrange [1e-2:1e-8]
set yrange [0:300]
set grid
set key inside left top box width -8 spacing 1.2 title "Calculated from \n {/CM-Typewriter vlt\\_HeatCap\\_P.c}"
set rmargin 4

plot 'vlt_HeatCap_P_00.050_DK0_1e-06.dat' using 2:9 with linespoints title " {/Times-Italic P} = {/Times-Italic P}_{sv} = 0.050 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.172 K) ", \
     'vlt_HeatCap_P_01.646_DK0_1e-06.dat' using 2:9 with linespoints title " {/Times-Italic P} = 1.644 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.157 K) ", \
     'vlt_HeatCap_P_07.328_DK0_1e-06.dat' using 2:9 with linespoints title " {/Times-Italic P} = 7.328 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.095 K) ", \
     'vlt_HeatCap_P_15.031_DK0_1e-06.dat' using 2:9 with linespoints title " {/Times-Italic P} = 15.031 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.998 K) ", \
     'vlt_HeatCap_P_18.180_DK0_1e-06.dat' using 2:9 with linespoints title " {/Times-Italic P} = 18.180 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.954 K) ", \
     'vlt_HeatCap_P_22.533_DK0_1e-06.dat' using 2:9 with linespoints title " {/Times-Italic P} = 22.533 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.889 K) ", \
     'vlt_HeatCap_P_25.868_DK0_1e-06.dat' using 2:9 with linespoints title " {/Times-Italic P} = 25.868 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.836 K) ", \
     -18.9102660394079*log(x) - 15.8225545948491 notitle, \
     -19.024807294548*log(x)  - 16.4107999003951 notitle, \
     -18.9299275918444*log(x) - 17.1139468917192 notitle, \
     -17.9293270703031*log(x) - 17.085805740948  notitle, \
     -17.3614842474275*log(x) - 17.2063164238303 notitle, \
     -16.5104564891659*log(x) - 17.2190939572875 notitle, \
     -15.8282446947989*log(x) - 16.5178636413462 notitle
