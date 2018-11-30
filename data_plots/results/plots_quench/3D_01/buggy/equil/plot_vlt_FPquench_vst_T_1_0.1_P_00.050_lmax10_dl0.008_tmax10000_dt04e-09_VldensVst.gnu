#filename: plot_vlt_FPquench_vst_T_1_0.1_P_00.050_lmax10_dl0.008_tmax10000_dt04e-09_VldensVst.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_vlt_FPquench_vst_T_1_0.1_P_00.050_lmax10_dl0.008_tmax10000_dt04e-09_VldensVst.ps'
#set size 1.0,1.1

set title "^4He Vortex-Loop Density vs Time {/Times-Italic t}/{/Times-Italic t}_0"

set ylabel "Vortex-Loop Density (units...)"
set xlabel "{/Times-Italic t}/{/Times-Italic t}_0 (unitless)"
#set xrange [0:4e-5]
#set xrange [0:4.38e-4]
#set yrange [0:3e-5]
set logscale xy
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot 'vlt_FPquench_vst_T_1_0.1_P_00.050_lmax10_dl0.008_tmax10000_dt04e-09.dat' using 2:9 with linespoints title "0.85 to 0.10"

#     '2DKT_FPquench_vst_T_1_0.1_lmax10_dl0.008_tmax10000_dt01e-05.dat' using 2:9 with linespoints title "1.00 to 0.10", \
#     '2DKT_FPquench_vst_T_0.98_0.1_lmax10_dl0.008_tmax10000_dt01e-05.dat' using 2:9 with linespoints title "0.98 to 0.10", \
#     '2DKT_FPquench_vst_lmax10_dl0.008_tmax10000_dt01e-05.dat' using 2:9 with linespoints title "0.95 to 0.10", \
#     '2DKT_FPquench_vst_T_0.95_0.1_lmax10_dl0.008_tmax10000_dt01e-05.dat' using 2:9 with linespoints title "0.95 to 0.10", \
#     '2DKT_FPquench_vst_T_0.9_0.1_lmax10_dl0.008_tmax10000_dt01e-05.dat' using 2:9 with linespoints title "0.90 to 0.10", \
#     '2DKT_FPquench_vst_T_0.85_0.1_lmax10_dl0.008_tmax10000_dt01e-05.dat' using 2:9 with linespoints title "0.85 to 0.10", \
#     '2DKT_FPquench_vst_T_0.8_0.1_lmax10_dl0.008_tmax10000_dt01e-05.dat' using 2:9 with linespoints title "0.80 to 0.10"
