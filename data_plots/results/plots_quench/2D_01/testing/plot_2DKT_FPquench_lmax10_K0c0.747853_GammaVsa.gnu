#filename: plot_2DKT_FPquench_lmax10_K0c0.747853_GammaVsa.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_2DKT_FPquench_lmax10_K0c0.747853_GammaVsa.ps'
#set size 1.0,1.1

set title "^4He Vortex-Pair Probability Density {/Symbol-Oblique G} vs Pair Separation {/Times-Italic a}/{/Times-Italic a}_0"

set ylabel "{/Symbol-Oblique G} (units...)"
set xlabel "{/Times-Italic a}/{/Times-Italic a}_0 (unitless)"
#set xrange [0:4e-5]
#set xrange [0:4.38e-4]
#set yrange [0:3e-5]
set logscale xy
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot '2DKT_FPquench_lmax10_dl0.008_tmax20600_dt00.00017.dat' using (exp($3)):6 with linespoints, \
     '2DKT_FPquench_lmax10_K0c0.747853.dat' using (exp($3)):6 with linespoints
#     'data.dat' using (exp($3)):5 with linespoints
