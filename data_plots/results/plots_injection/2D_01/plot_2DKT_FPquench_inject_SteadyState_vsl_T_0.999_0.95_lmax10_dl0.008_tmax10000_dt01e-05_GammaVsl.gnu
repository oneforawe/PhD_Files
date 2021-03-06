#filename: plot_2DKT_FPquench_inject_SteadyState_vsl_T_0.999_0.95_lmax10_dl0.008_tmax10000_dt01e-05_GammaVsl.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_2DKT_FPquench_inject_SteadyState_vsl_T_0.999_0.95_lmax10_dl0.008_tmax10000_dt01e-05_GammaVsl.ps'
#set size 1.0,1.1

set title "^4He Vortex-Pair Distribution {/Symbol-Oblique G} vs Pair Separation Parameter {/Times-Italic l}"

set ylabel "{/Symbol-Oblique G} (units...)"
set xlabel "{/Times-Italic l} (unitless)"
set xrange [0:4]
set yrange [1e-15:1]
set logscale y
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot '2DKT_FPquench_inject_SteadyState_vsl_T_1_lmax10_dl0.008_linj3_InjConst1_dGdl00.dat' using 3:6 with linespoints title " At Tc (w/ injection) ", \
     '2DKT_FPquench_vsl_T_1_0.95_lmax10_dl0.008_tmax10000_dt01e-05_equil.dat' using 3:6 with linespoints title " At Tc (no injection) ", \
     '2DKT_FPquench_inject_SteadyState_vsl_T_1_lmax10_dl0.008_linj3_InjConst1_dGdl00.01.dat' using 3:6 with linespoints title " 1.0 of Tc (w/ injection) ", \
     '2DKT_FPquench_inject_SteadyState_vsl_T_1_lmax10_dl0.008_linj3_InjConst1_dGdl0-0.001.dat' using 3:6 with linespoints title " 1.0 of Tc (w/ injection) "
#     '2DKT_FPquench_inject_SteadyState_vsl_T_0.75_lmax10_dl0.008_linj3_InjConst1_dGdl00.dat' using 3:6 with linespoints title " 0.75 of Tc (w/ injection) ", \
#     '2DKT_FPquench_vsl_T_0.75_0.95_lmax10_dl0.008_tmax10000_dt01e-05_equil.dat' using 3:6 with linespoints title " 0.75 of Tc (no injection) "
