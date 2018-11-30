#filename: plot_vlt_FPquench_vsl_T_1_0.1_P_00.050_lmax10_dl0.008_tmax10000_dt04e-09_GammaVsa.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_vlt_FPquench_vsl_T_1_0.1_P_00.050_lmax10_dl0.008_tmax10000_dt04e-09_GammaVsa.ps'
#set size 1.0,1.1

set title "^4He Vortex-Loop Distribution {/Symbol-Oblique G} vs Loop Diameter {/Times-Italic a}/{/Times-Italic a}_0"

set ylabel "{/Symbol-Oblique G} (units...)"
set xlabel "{/Times-Italic a}/{/Times-Italic a}_0 (unitless)"
#set xrange [1:1.1]
#set yrange [0:3e-5]
set logscale y
#set logscale xy
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

#plot 'vlt_FPquench_vsl_T_1_0.1_P_00.050_lmax10_dl0.008_tmax10000_dt04e-09.dat' using (exp($3)):7 with linespoints
plot 'vlt_FPquench_vsl_T_1_0.1_P_00.050_lmax10_dl0.008_tmax10000_dt04e-09.dat' using 3:7 with linespoints, \
     'vlt_FPquench_vsl_T_1_0.1_P_00.050_lmax10_dl0.008_tmax10000_dt01e-12.dat' using 3:7 with linespoints
