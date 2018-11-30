#filename: plot_HeatCapVsTempV.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_HeatCapVsTempV.ps'

set title "^4He Molar Specific Heat Capacity at Constant Pressure {/Times-Italic c_P} \n vs Temperature Variable {/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l}"

set ylabel "{/Times-Italic c_P} (J mol^{-1} K^{-1})"
set xlabel "{/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l} (unitless)"
set xrange [1:10**(-8)]
set yrange [0:150]
set logscale x
set grid
set key inside left top box width -6.5 height 0.5 spacing 1.1 title " Calculated using \n {/CM-Typewriter vlt\\_HeatCap.c} "
set rmargin 2

plot 'vlt_HeatCap_DK0_1e-8.dat' using 1:2 title ' DK0 = 1e-8 ' with linespoints, \
     'vlt_HeatCap_DK0_1e-7.dat' using 1:2 title ' DK0 = 1e-7 ' with linespoints, \
     'vlt_HeatCap_DK0_1e-6.dat' using 1:2 title ' DK0 = 1e-6 ' with linespoints, \
     'vlt_HeatCap_DK0_1e-5.dat' using 1:2 title ' DK0 = 1e-5 ' with linespoints, \
     'vlt_HeatCap_DK0_1e-4.dat' using 1:2 title ' DK0 = 1e-4 ' with linespoints, \
     'vlt_HeatCap_DK0_1e-3.dat' using 1:2 title ' DK0 = 1e-3 ' with linespoints, \
     'vlt_HeatCap_DK0_1e-2.dat' using 1:2 title ' DK0 = 1e-2 ' with linespoints, \
     'vlt_HeatCap_DK0_1e-1.dat' using 1:2 title ' DK0 = 1e-1 ' with linespoints, \
     -13.7720820955*log10(x)-5.67096106256 title ' fit = -14*log10({/Symbol-Oblique t})-6', \
     "<echo 3.65982e-07  82.9736;  echo 0.0080613  23.1622" title ' fit points'
