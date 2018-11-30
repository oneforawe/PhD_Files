#filename: plot_2Dvpt_constT_lmax100_dl0.01_Op1_SigsSigVsl.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_2Dvpt_constT_lmax100_dl0.01_Op1_SigsSigVsl.ps'
#set size 1.0,1.1

set title "^4He Superfluid Fraction {/Symbol-Oblique s}_s/{/Symbol-Oblique s} vs Pair Separation Parameter {/Times-Italic l}"

set ylabel "{/Symbol-Oblique s}_s/{/Symbol-Oblique s} (unitless)"
set xlabel "{/Times-Italic l} (unitless)"
#set xrange [0:1]
#set yrange [5:8]
#set logscale y
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot '2Dvpt_constT_0.205672_lmax100_dl0.01_Op1.dat' using 3:5 with linespoints title " 0.206*Tc ", \
     '2Dvpt_constT_0.9_lmax100_dl0.01_Op1.dat'      using 3:5 with linespoints title " 0.900*Tc ", \
     '2Dvpt_constT_0.99_lmax100_dl0.01_Op1.dat'     using 3:5 with linespoints title " 0.990*Tc ", \
     '2Dvpt_constT_0.999_lmax100_dl0.01_Op1.dat'    using 3:5 with linespoints title " 0.999*Tc ", \
     '2Dvpt_constT_1.00126_lmax100_dl0.01_Op1.dat'  using 3:5 with linespoints title " 1.001*Tc ", \
     '2Dvpt_constT_1.01_lmax100_dl0.01_Op1.dat'     using 3:5 with linespoints title " 1.010*Tc ", \
     '2Dvpt_constT_1.1_lmax100_dl0.01_Op1.dat'      using 3:5 with linespoints title " 1.100*Tc ", \
     '2Dvpt_constT_2_lmax100_dl0.01_Op1.dat'        using 3:5 with linespoints title " 2.000*Tc "
