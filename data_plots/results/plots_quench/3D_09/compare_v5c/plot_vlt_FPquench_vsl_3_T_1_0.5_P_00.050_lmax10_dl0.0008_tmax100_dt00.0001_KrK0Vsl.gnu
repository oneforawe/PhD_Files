#filename: plot_vlt_FPquench_vsl_3_T_1_0.5_P_00.050_lmax10_dl0.0008_tmax100_dt00.0001_KrK0Vsl.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_vlt_FPquench_vsl_3_T_1_0.5_P_00.050_lmax10_dl0.0008_tmax100_dt00.0001_KrK0Vsl.ps'
#set size 1.0,1.1

set title "^4He Superfluid Fraction {/Times-Italic K}_r/{/Times-Italic K}_0 vs Loop Diameter Scale {/Times-Italic l}"

set ylabel "{/Times-Italic K}_r/{/Times-Italic K}_0 (unitless)"
set xlabel "{/Times-Italic l} (unitless)"
#set xrange [0:5]
#set yrange [1e-20:1]
#set xrange [0:3]
#set yrange [1e-9:1]
#set xrange [0:0.25]
#set yrange [0.01:0.1]
#set xrange [0:2]
#set yrange [1e-10:1]
#set yrange [1e-10:1]
#set xrange [0:4e-5]
#set xrange [0:4.38e-4]
#set yrange [1e-5:1]
set logscale y
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot 'vlt_FPquench_vsl_3_T_1_0.5_P_00.050_lmax10_dl0.0008_tmax100_dt00.0001.dat' using 3:6 with linespoints
