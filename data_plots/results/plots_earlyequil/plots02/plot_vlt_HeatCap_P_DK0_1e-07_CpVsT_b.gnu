#filename: plot_vlt_HeatCap_P_DK0_1e-07_CpVsT_b.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_vlt_HeatCap_P_DK0_1e-07_CpVsT_b.ps'

set title "^4He Molar Specific Heat Capacity at Constant Pressure {/Times-Italic c_P} \n vs Exponentiated Temperature Variable {/Symbol-Oblique t}^{-{/Symbol-Oblique a}}"

set ylabel "{/Times-Italic c_P} (J mol^{-1} K^{-1})"
set xlabel "{/Symbol-Oblique t}^{-{/Symbol-Oblique a}} = [{/Times-Italic 1 - T / T}_{/Symbol-Oblique l}{/Times-Italic (P)}]^{0.0150649} (unitless)"
#set logscale x
set xrange [0.9998:1.00005]
set yrange [100:400]
set grid
set key inside left top box width -4 spacing 1.2 title "Calculated from \n {/CM-Typewriter vlt\\_HeatCap\\_P.c}"
set rmargin 2

set sample 10000


plot 'vlt_HeatCap_P_00.050_DK0_1e-07.dat' using (1-$1):2 with points title " {/Times-Italic P} = {/Times-Italic P}_{sv} = 0.050 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.172 K) ", \
     'vlt_HeatCap_P_01.646_DK0_1e-07.dat' using (1-$1):2 with points title " {/Times-Italic P} = 1.644 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.157 K) ", \
     'vlt_HeatCap_P_07.328_DK0_1e-07.dat' using (1-$1):2 with points title " {/Times-Italic P} = 7.328 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.095 K) ", \
     'vlt_HeatCap_P_15.031_DK0_1e-07.dat' using (1-$1):2 with points title " {/Times-Italic P} = 15.031 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.998 K) ", \
     'vlt_HeatCap_P_18.180_DK0_1e-07.dat' using (1-$1):2 with points title " {/Times-Italic P} = 18.180 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.954 K) ", \
     'vlt_HeatCap_P_22.533_DK0_1e-07.dat' using (1-$1):2 with points title " {/Times-Italic P} = 22.533 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.889 K) ", \
     'vlt_HeatCap_P_25.868_DK0_1e-07.dat' using (1-$1):2 with points title " {/Times-Italic P} = 25.868 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.836 K) ", \
     -1421.41776427837*(1-x)**0.0150649 + 1396.0719948076 notitle, \
     -1429.87948666493*(1-x)**0.0150649 + 1403.90623666182 notitle, \
     -1425.00044310711*(1-x)**0.0150649 + 1398.14824914441 notitle, \
     -1354.17794882653*(1-x)**0.0150649 + 1327.24949408352 notitle, \
     -1309.4203735063*(1-x)**0.0150649  + 1283.22903912049 notitle, \
     -1238.3514407305*(1-x)**0.0150649  + 1213.31118757252 notitle, \
     -1196.33128220878*(1-x)**0.0150649 + 1171.07642207644 notitle, \
     -42.9309316155699*log10(1-x) - 13.9036180158113 notitle, \
     -43.1884331517886*log10(1-x) - 14.4707240110324 notitle, \
     -43.0415626666153*log10(1-x) - 15.3909086004666 notitle, \
     -40.9069885145288*log10(1-x) - 16.0551780694881 notitle, \
     -39.5489034994062*log10(1-x) - 15.6532134023518 notitle, \
     -37.4023206551689*log10(1-x) - 15.0738422848683 notitle, \
     -36.1353190280632*log10(1-x) - 15.6352286131365 notitle
