#filename: plot_vlt_FPquench_vsl_3_P_00.050_lmax10_dl0.008_equil_GammaVsa.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_vlt_FPquench_vsl_3_P_00.050_lmax10_dl0.008_equil_GammaVsa.ps'
#set size 1.0,1.1

set title "Equilibrium ^4He Superfluid Fraction {/Symbol-Oblique r}/{/Symbol-Oblique r}_s vs Loop Diameter Scale {/Times-Italic l}"

set ylabel "{/Symbol-Oblique r}/{/Symbol-Oblique r}_s (unitless)"
set xlabel "{/Times-Italic l} (unitless)"
#set xrange [0:4e-5]
#set xrange [0:4.38e-4]
set yrange [1e-5:1]
set logscale y
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot 'vlt_FPquench_vsl_3_T_1_0.1_P_00.050_lmax10_dl0.008_tmax0.01_dt00.0001_equil.dat' using 3:6 with linespoints title " {/Times-Italic T}/{/Times-Italic T}_c = 1.00 ", \
     'vlt_FPquench_vsl_3_T_0.98_0.1_P_00.050_lmax10_dl0.008_tmax0.01_dt00.0001_equil.dat' using 3:6 with linespoints title " {/Times-Italic T}/{/Times-Italic T}_c = 0.98 ", \
     'vlt_FPquench_vsl_3_T_0.95_0.1_P_00.050_lmax10_dl0.008_tmax0.01_dt00.0001_equil.dat' using 3:6 with linespoints title " {/Times-Italic T}/{/Times-Italic T}_c = 0.95 ", \
     'vlt_FPquench_vsl_3_T_0.75_0.1_P_00.050_lmax10_dl0.008_tmax0.01_dt00.0001_equil.dat' using 3:6 with linespoints title " {/Times-Italic T}/{/Times-Italic T}_c = 0.75 ", \
     'vlt_FPquench_vsl_3_T_0.5_0.1_P_00.050_lmax10_dl0.008_tmax0.01_dt00.0001_equil.dat' using 3:6 with linespoints title " {/Times-Italic T}/{/Times-Italic T}_c = 0.50 "

#     'vlt_FPquench_vsl_3_T_1_0.1_P_00.050_lmax10_dl0.008_tmax10000_dt00.0001_equil.dat' using 3:5 with linespoints title " K/K0 {/Times-Italic T}/{/Times-Italic T}_c = 1.00 ", \
#     'vlt_FPquench_vsl_3_T_0.95_0.1_P_00.050_lmax10_dl0.008_tmax10000_dt00.0001_equil.dat' using 3:5 with linespoints title " K/K0 {/Times-Italic T}/{/Times-Italic T}_c = 0.95 ", \
#     'vlt_FPquench_vsl_3_T_0.75_0.1_P_00.050_lmax10_dl0.005_tmax10000_dt00.0001_equil.dat' using 3:6 with linespoints title " {/Times-Italic T}/{/Times-Italic T}_c = 0.75 ", \
#     'vlt_FPquench_vsl_3_T_0.75_0.1_P_00.050_lmax10_dl0.008_tmax10000_dt00.0001_equil.dat' using 3:5 with linespoints title " K/K0 {/Times-Italic T}/{/Times-Italic T}_c = 0.75 ", \
#     'vlt_FPquench_vsl_3_T_0.75_0.1_P_00.050_lmax10_dl0.005_tmax10000_dt00.0001_equil.dat' using 3:5 with linespoints title " K/K0 {/Times-Italic T}/{/Times-Italic T}_c = 0.75 ", \
#     'vlt_FPquench_vsl_3_T_0.5_0.1_P_00.050_lmax10_dl0.008_tmax10000_dt00.0001_equil.dat' using 3:5 with linespoints title " K/K0 {/Times-Italic T}/{/Times-Italic T}_c = 0.50 "
