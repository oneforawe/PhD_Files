#filename: plot_3Dvlt_quench_vsl_T_1_0.7_P_00.050_lmax10_dl0.001_tmax1_dt5e-07_BC1_Kvsl.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_3Dvlt_quench_vsl_T_1_0.7_P_00.050_lmax10_dl0.001_tmax1_dt5e-07_BC1_Kvsl.ps'
#set size 1.0,1.1

set title "^4He Coupling Parameter {/Times-Italic K} vs Length Scale {/Times-Italic l}"

set ylabel "{/Times-Italic K} (unitless)"
set xlabel "{/Times-Italic l} (unitless)"
#set yrange [0:1000]
#set xrange [0:4e-5]
#set xrange [0:4.38e-4]
#set yrange [1e-30:1]
#set logscale y
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot '3Dvlt_quench_vsl_T_1_0.7_P_00.050_lmax10_dl0.001_tmax1_dt5e-07_BC1.dat' using 3:4 with linespoints
