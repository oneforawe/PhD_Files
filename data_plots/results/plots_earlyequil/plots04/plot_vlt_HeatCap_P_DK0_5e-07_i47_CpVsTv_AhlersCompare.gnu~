#filename: plot_vlt_HeatCap_P_DK0_5e-07_i47_CpVsTv_AhlersCompare.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_vlt_HeatCap_P_DK0_5e-07_i47_CpVsTv_AhlersCompare.ps'

set title "^4He Molar Specific Heat Capacity at Constant Pressure {/Times-Italic c_P} \n vs Temperature Variable {/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l}"

set ylabel "{/Times-Italic c_P} (J mol^{-1} K^{-1})"
set xlabel "{/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l}{/Times-Italic (P)} (unitless)"
set logscale x
set xrange [1:1e-8]
#set xrange [1e-2:1e-8]
set mxtics 10
set yrange [0:300]
set grid
set key inside left top box width -8 spacing 1.2 title "Calculated from \n {/CM-Typewriter vlt\\_HeatCap\\_P.c}"
set rmargin 4

plot 'vlt_HeatCap_P_00.050_DK0_5e-07_i47.dat' using 2:9 with linespoints title " {/Times-Italic P} = {/Times-Italic P}_{sv} = 0.050 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.172 K) ", \
     'vlt_HeatCap_P_01.646_DK0_5e-07_i47.dat' using 2:9 with linespoints title " {/Times-Italic P} = 1.644 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.157 K) ", \
     'vlt_HeatCap_P_07.328_DK0_5e-07_i47.dat' using 2:9 with linespoints title " {/Times-Italic P} = 7.328 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.095 K) ", \
     'vlt_HeatCap_P_15.031_DK0_5e-07_i47.dat' using 2:9 with linespoints title " {/Times-Italic P} = 15.031 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.998 K) ", \
     'vlt_HeatCap_P_18.180_DK0_5e-07_i47.dat' using 2:9 with linespoints title " {/Times-Italic P} = 18.180 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.954 K) ", \
     'vlt_HeatCap_P_22.533_DK0_5e-07_i47.dat' using 2:9 with linespoints title " {/Times-Italic P} = 22.533 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.889 K) ", \
     'vlt_HeatCap_P_25.868_DK0_5e-07_i47.dat' using 2:9 with linespoints title " {/Times-Italic P} = 25.868 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.836 K) ", \
     -1423.94920626784*x**0.0150649  + 1398.33503103225 notitle, \
     -1432.89744603034*x**0.0150649  + 1406.59321021477 notitle, \
     -1427.01434591532*x**0.0150649  + 1399.89715998998 notitle, \
     -1352.37023755605*x**0.0150649  + 1325.71266994426 notitle, \
     -1310.16511033501*x**0.0150649  + 1283.59606938083 notitle, \
     -1246.29923838762*x**0.0150649  + 1220.13906155327 notitle, \
     -1196.63376643744*x**0.0150649  + 1171.30717773604 notitle
