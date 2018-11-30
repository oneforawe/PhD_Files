#filename: plot_2Dvpt_constT_lmax10_dl0.008_SigsSigVsl.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_2Dvpt_constT_lmax10_dl0.008_SigsSigVsl.ps'
#set size 1.0,1.1

set title "^4He Superfluid Fraction {/Symbol-Oblique s}_s/{/Symbol-Oblique s} vs Pair Separation Parameter {/Times-Italic l}"

set ylabel "{/Symbol-Oblique s}_s/{/Symbol-Oblique s} (unitless)"
set xlabel "{/Times-Italic l} (unitless)"
#set xrange [0:1]
set yrange [0:1]
#set logscale y
set grid
set key outside right box width -2 height 0.5 spacing 1.2 title " 2D calcs "
#set rmargin 35

plot '2Dvpt_constT_1_lmax10_dl0.008.dat'   using 3:5 with linespoints title " 1.0*Tc ", \
     '2Dvpt_constT_0.9_lmax10_dl0.008.dat' using 3:5 with linespoints title " 0.9*Tc ", \
     '2Dvpt_constT_0.8_lmax10_dl0.008.dat' using 3:5 with linespoints title " 0.8*Tc ", \
     '2Dvpt_constT_0.7_lmax10_dl0.008.dat' using 3:5 with linespoints title " 0.7*Tc ", \
     '2Dvpt_constT_0.6_lmax10_dl0.008.dat' using 3:5 with linespoints title " 0.6*Tc ", \
     '2Dvpt_constT_0.5_lmax10_dl0.008.dat' using 3:5 with linespoints title " 0.5*Tc ", \
     '2Dvpt_constT_0.4_lmax10_dl0.008.dat' using 3:5 with linespoints title " 0.4*Tc ", \
     '2Dvpt_constT_0.3_lmax10_dl0.008.dat' using 3:5 with linespoints title " 0.3*Tc ", \
     '2Dvpt_constT_0.2_lmax10_dl0.008.dat' using 3:5 with linespoints title " 0.2*Tc ", \
     '2Dvpt_constT_0.1_lmax10_dl0.008.dat' using 3:5 with linespoints title " 0.1*Tc "
