#filename: plot_2DKT_FPquench_vsl_lmax10_dl0.008_tmax10000_dt01e-05_GammaVsa.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_2DKT_FPquench_vsl_lmax10_dl0.008_tmax10000_dt01e-05_GammaVsa.ps'
#set size 1.0,1.1

set title "^4He Vortex-Pair Distribution {/Symbol-Oblique G} vs Pair Separation Parameter {/Times-Italic l}"

set ylabel "{/Symbol-Oblique G} (units...)"
set xlabel "{/Times-Italic l} (unitless)"
# Look at the weird plot from 2Dquench_ProbDensG.dat (2Dquench.c)
#set xrange [0:0.5]
#set yrange [1e-4:1e-3]
# wider look
#set xrange [0:1]
#set yrange [1e-6:1e-3]
# yet wider
set xrange [0:2]
set yrange [1e-8:1e-3]
# different look
#set xrange [2:2.5]
#set yrange [1e-9:1e-7]
# different look
#set xrange [0.5:1.5]
#set yrange [1e-7:1e-4]
set logscale y
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot '2DKT_FPquench_vsl_T_0.95_0.85_lmax10_dl0.008_tmax10000_dt01e-05.dat' using 3:6 with linespoints title " 0.95 to 0.85 of Tc ", \
     '2DKT_FPquench_vsl_T_0.999_0.95_lmax10_dl0.008_tmax10000_dt01e-05.dat' using 3:6 with linespoints title " 0.999 to 0.95 of Tc "
