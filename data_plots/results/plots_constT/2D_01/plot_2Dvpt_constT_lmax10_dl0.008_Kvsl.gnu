#filename: plot_2Dvpt_constT_lmax10_dl0.008_Kvsl.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_2Dvpt_constT_lmax10_dl0.008_Kvsl.ps'
#set size 1.0,1.1

set title "^4He Superfluid-Temperature Ratio {/Times-Italic K} vs Pair Separation Parameter {/Times-Italic l}"

set ylabel "{/Times-Italic K} (unitless)"
set xlabel "{/Times-Italic l} (unitless)"
#set xrange [0:1]
#set yrange [5:8]
#set logscale y
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot '2Dvpt_constT_1_lmax10_dl0.008.dat'   using 3:4 with linespoints title " 1.0*Tc ", \
     '2Dvpt_constT_0.9_lmax10_dl0.008.dat' using 3:4 with linespoints title " 0.9*Tc ", \
     '2Dvpt_constT_0.8_lmax10_dl0.008.dat' using 3:4 with linespoints title " 0.8*Tc ", \
     '2Dvpt_constT_0.7_lmax10_dl0.008.dat' using 3:4 with linespoints title " 0.7*Tc ", \
     '2Dvpt_constT_0.6_lmax10_dl0.008.dat' using 3:4 with linespoints title " 0.6*Tc ", \
     '2Dvpt_constT_0.5_lmax10_dl0.008.dat' using 3:4 with linespoints title " 0.5*Tc ", \
     '2Dvpt_constT_0.4_lmax10_dl0.008.dat' using 3:4 with linespoints title " 0.4*Tc ", \
     '2Dvpt_constT_0.3_lmax10_dl0.008.dat' using 3:4 with linespoints title " 0.3*Tc ", \
     '2Dvpt_constT_0.2_lmax10_dl0.008.dat' using 3:4 with linespoints title " 0.2*Tc ", \
     '2Dvpt_constT_0.1_lmax10_dl0.008.dat' using 3:4 with linespoints title " 0.1*Tc "
