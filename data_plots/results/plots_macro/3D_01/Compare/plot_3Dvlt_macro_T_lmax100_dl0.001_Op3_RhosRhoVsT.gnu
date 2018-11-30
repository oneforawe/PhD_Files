#filename: plot_3Dvlt_macro_T_lmax100_dl0.001_Op3_RhosRhoVsT.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_3Dvlt_macro_T_lmax100_dl0.001_Op3_RhosRhoVsT.ps'
#set size 1.0,1.1

set title "^4He Superfluid Fraction {/Symbol-Oblique r}_s/{/Symbol-Oblique r} vs Temperature Fraction {/Times-Italic T}/{/Times-Italic T}_c"

set ylabel "{/Symbol-Oblique r}_s/{/Symbol-Oblique r} (unitless)"
set xlabel "{/Times-Italic T}/{/Times-Italic T}_c (unitless)"
set xrange [0.2:1.4]
#set xrange [0:1]
#set yrange [5:8]
#set logscale y
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot '3Dvlt_macro_T_Cc1.20_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 1.20 ", \
     '3Dvlt_macro_T_Cc1.10_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 1.10 ", \
     '3Dvlt_macro_T_Cc1.06_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 1.06 ", \
     '3Dvlt_macro_T_Cc1.05_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 1.05 ", \
     '3Dvlt_macro_T_Cc1.04_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 1.04 ", \
     '3Dvlt_macro_T_Cc1.03_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 1.03 ", \
     '3Dvlt_macro_T_Cc1.02_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 1.02 ", \
     '3Dvlt_macro_T_Cc1.01_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 1.01 ", \
     '3Dvlt_macro_T_Cc1.00_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 1.00 ", \
     '3Dvlt_macro_T_Cc0.99_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.99 ", \
     '3Dvlt_macro_T_Cc0.98_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.98 ", \
     '3Dvlt_macro_T_Cc0.97_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.97 ", \
     '3Dvlt_macro_T_Cc0.90_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.90 ", \
     '3Dvlt_macro_T_Cc0.80_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.80 ", \
     '3Dvlt_macro_T_Cc0.70_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.70 ", \
     '3Dvlt_macro_T_Cc0.60_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.60 ", \
     '3Dvlt_macro_T_Cc0.55_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.55 ", \
     '3Dvlt_macro_T_Cc0.50_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.50 ", \
     '3Dvlt_macro_T_Cc0.40_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.40 ", \
     '3Dvlt_macro_T_Cc0.30_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.30 ", \
     '3Dvlt_macro_T_Cc0.20_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.20 ", \
     '3Dvlt_macro_T_Cc0.10_lmax100_dl0.001_Op3.dat' using 1:6 with linespoints title " Cc = 0.10 "
