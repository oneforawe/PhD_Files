#filename: plot_3Dvlt_constT_Cc1.11_lmax10_dl0.001_Op1_GammaVsl.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_3Dvlt_constT_Cc1.11_lmax10_dl0.001_Op1_GammaVsl.ps'
#set size 1.0,1.1

set title "^4He Vortex-Pair Distribution {/Symbol-Oblique G} vs Pair Separation Parameter {/Times-Italic l}"

set ylabel "{/Symbol-Oblique G} (units...)"
set xlabel "{/Times-Italic l} (unitless)"
#set xrange [0:0.005]
#set xrange [0:0.02]
#set yrange [1e-20:1]
set logscale y
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot '3Dvlt_constT_1_Cc1.11_lmax10_dl0.001_Op1.dat'   using 3:7 with linespoints title " 1.0*Tc ", \
     '3Dvlt_constT_0.9_Cc1.11_lmax10_dl0.001_Op1.dat' using 3:7 with linespoints title " 0.9*Tc ", \
     '3Dvlt_constT_0.8_Cc1.11_lmax10_dl0.001_Op1.dat' using 3:7 with linespoints title " 0.8*Tc ", \
     '3Dvlt_constT_0.7_Cc1.11_lmax10_dl0.001_Op1.dat' using 3:7 with linespoints title " 0.7*Tc ", \
     '3Dvlt_constT_0.6_Cc1.11_lmax10_dl0.001_Op1.dat' using 3:7 with linespoints title " 0.6*Tc ", \
     '3Dvlt_constT_0.5_Cc1.11_lmax10_dl0.001_Op1.dat' using 3:7 with linespoints title " 0.5*Tc ", \
     '3Dvlt_constT_0.4_Cc1.11_lmax10_dl0.001_Op1.dat' using 3:7 with linespoints title " 0.4*Tc ", \
     '3Dvlt_constT_0.3_Cc1.11_lmax10_dl0.001_Op1.dat' using 3:7 with linespoints title " 0.3*Tc ", \
     '3Dvlt_constT_0.2_Cc1.11_lmax10_dl0.001_Op1.dat' using 3:7 with linespoints title " 0.2*Tc ", \
     '3Dvlt_constT_0.1_Cc1.11_lmax10_dl0.001_Op1.dat' using 3:7 with linespoints title " 0.1*Tc "
