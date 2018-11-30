#filename: plot_3Dvlt_macro_T_lmax100_dl0.001_Op3_RhosRhoVsT.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_3Dvlt_macro_T_lmax100_dl0.001_Op3_RhosRhoVsT.ps'
#set size 1.0,1.1

set title "^4He Superfluid Fraction {/Symbol-Oblique r}_s/{/Symbol-Oblique r} vs Temperature Fraction {/Times-Italic T}/{/Times-Italic T}_c"

set ylabel "{/Symbol-Oblique r}_s/{/Symbol-Oblique r} (unitless)"
set xlabel "{/Times-Italic T}/{/Times-Italic T}_c (unitless)"
#set xrange [0:1]
#set yrange [5:8]
#set logscale y
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot '3Dvlt_macro_T_Cc0.55_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.55 ", \
     '3Dvlt_macro_T_Cc0.50_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.50 ", \
     '3Dvlt_macro_T_Cc0.45_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.45 ", \
     '3Dvlt_macro_T_Cc0.40_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.40 ", \
     '3Dvlt_macro_T_Cc0.35_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.35 ", \
     '3Dvlt_macro_T_Cc0.30_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.30 ", \
     '3Dvlt_macro_T_Cc0.25_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.25 ", \
     '3Dvlt_macro_T_Cc0.20_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.20 ", \
     '3Dvlt_macro_T_Cc0.15_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.15 ", \
     '3Dvlt_macro_T_Cc0.10_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.10 ", \
     '3Dvlt_macro_T_Cc0.05_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.05 "
