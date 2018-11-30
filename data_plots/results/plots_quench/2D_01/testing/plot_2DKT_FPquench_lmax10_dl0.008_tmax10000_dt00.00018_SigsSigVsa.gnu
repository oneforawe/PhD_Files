#filename: plot_2DKT_FPquench_lmax10_dl0.008_tmax10000_dt00.00018_SigsSigVsa.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_2DKT_FPquench_lmax10_dl0.008_tmax10000_dt00.00018_SigsSigVsa.ps'
#set size 1.0,1.1

set title "^4He Superfluid Fraction {/Symbol-Oblique s}_s/{/Symbol-Oblique s} vs Pair Separation {/Times-Italic a}/{/Times-Italic a}_0"

set ylabel "{/Symbol-Oblique s}_s/{/Symbol-Oblique s} (unitless)"
set xlabel "{/Times-Italic a}/{/Times-Italic a}_0 (unitless)"
#set xrange [0:4e-5]
#set xrange [0:4.38e-4]
#set yrange [0:3e-5]
set logscale x
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot '2DKT_FPquench_glitch_lmax10_dl0.008_tmax10000_dt00.00018.dat' using (exp($3)):5 with linespoints lc rgb "green" title "Weird eqn, glitch at start", \
     '2DKT_FPquench_lmax10_dl0.008_tmax10000_dt00.00018.dat' using (exp($3)):5 with linespoints lc rgb "red" title "Weird eqn, corrected glitch", \
     '2DKT_FPquench_lmax10_dl0.008_tmax10000_dt01e-05.dat' using (exp($3)):5 with linespoints pt 4 title "Correct equation"
