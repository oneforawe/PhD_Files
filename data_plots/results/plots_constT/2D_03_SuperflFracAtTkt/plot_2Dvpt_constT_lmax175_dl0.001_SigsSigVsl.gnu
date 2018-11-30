#filename: plot_2Dvpt_constT_lmax175_dl0.001_SigsSigVsl.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_2Dvpt_constT_lmax175_dl0.001_SigsSigVsl.ps'
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

plot '2Dvpt_constT_1_lmax175_dl0.001_Op1.dat' using 3:5 with linespoints title " 1.0*Tc ", \
     (2.0/pi)/0.747853

# K*/K0 = (2/pi)/0.747853