#filename: plot_3Dvlt_macro_T_Cc1.11_lmax100_dl0.001_Op3_GammaVsT.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_3Dvlt_macro_T_Cc1.11_lmax100_dl0.001_Op3_GammaVsT.ps'
#set size 1.0,1.1

set title "^4He Vortex-Loop Distribution {/Symbol-Oblique G} vs Temperature Fraction {/Times-Italic T}/{/Times-Italic T}_c"

set ylabel "{/Symbol-Oblique G} (units...)"
set xlabel "{/Times-Italic T}/{/Times-Italic T}_c (unitless)"
#set xrange [0:0.005]
#set xrange [0:0.02]
#set yrange [1e-20:1]
set logscale y
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot '3Dvlt_macro_T_Cc1.11_lmax1_dl0.001_Op3.dat'   using 1:7 with linespoints title " lmax = 1 ", \
     '3Dvlt_macro_T_Cc1.11_lmax10_dl0.001_Op3.dat'  using 1:7 with linespoints title " lmax = 10 ", \
     '3Dvlt_macro_T_Cc1.11_lmax100_dl0.001_Op3.dat' using 1:7 with linespoints title " lmax = 100 "
