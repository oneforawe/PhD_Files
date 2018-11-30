#filename: plot_3Dvlt_macro_cp_P_Dk05e-07_i47_RhosoRhoVsT.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_3Dvlt_macro_cp_P_Dk05e-07_i47_RhosoRhoVsT.ps'
set size 0.85,0.85

set title "^4He Superfluid Fraction {/Symbol-Oblique r}_s / {/Symbol-Oblique r} vs Temperature Fraction {/Times-Italic T} / {/Times-Italic T}_{/Symbol-Oblique l} "

set ylabel "{/Symbol-Oblique r}_s / {/Symbol-Oblique r} (unitless)"
set xlabel "{/Times-Italic T} / {/Times-Italic T}_{/Symbol-Oblique l}({/Times-Italic P}) (unitless)"
set xrange [0:1]
#set xrange [0.95:1]
#set xrange [0.999:1]
#set xrange [0.99984:1]

## Version 1 ##
#set xrange [0.9998:1]
#set xtics 0.00005

## Version 2 ##
#set xrange [0.9999:0.99998]
#set yrange [0.001:0.0055]
#set xtics 0.00005

#set grid
set key outside right box width -7 height 0.5 spacing 1.2 title " Vortex Loop Theory vs Data \n of Greywall and Ahlers \n\n {/Times-Italic P}             ({/Times-Italic T}_{/Symbol-Oblique l})           "
set rmargin 25

#set xtics 0.99984 0.00004 1.0
#set mxtics 2
set pointsize 0.85

plot '3Dvlt_macro_cp_P00.050_DK05e-07_i47.dat' using 1:9 with lines title "  0.05 bar  (2.172 K) ", \
     '3Dvlt_macro_cp_P07.270_DK05e-07_i47.dat' using 1:9 with lines title "  7.27 bar  (2.096 K) ", \
     '3Dvlt_macro_cp_P12.130_DK05e-07_i47.dat' using 1:9 with lines title " 12.13 bar  (2.036 K) ", \
     '3Dvlt_macro_cp_P18.060_DK05e-07_i47.dat' using 1:9 with lines title " 18.06 bar  (1.956 K) ", \
     '3Dvlt_macro_cp_P24.100_DK05e-07_i47.dat' using 1:9 with lines title " 24.10 bar  (1.865 K) ", \
     '3Dvlt_macro_cp_P29.090_DK05e-07_i47.dat' using 1:9 with lines title " 29.09 bar  (1.782 K) ", \
     'GreywallAhlers_a_checked.dat' using (1-$1):3 with linespoints linestyle 1 lc rgb "black" title " Data ", \
     'GreywallAhlers_a_checked.dat' using (1-$4):6 with linespoints linestyle 1 lc rgb "black" notitle, \
     'GreywallAhlers_a_checked.dat' using (1-$7):9 with linespoints linestyle 1 lc rgb "black" notitle, \
     'GreywallAhlers_b_checked.dat' using (1-$1):3 with linespoints linestyle 1 lc rgb "black" notitle, \
     'GreywallAhlers_b_checked.dat' using (1-$4):6 with linespoints linestyle 1 lc rgb "black" notitle, \
     'GreywallAhlers_b_checked.dat' using (1-$7):9 with linespoints linestyle 1 lc rgb "black" notitle
