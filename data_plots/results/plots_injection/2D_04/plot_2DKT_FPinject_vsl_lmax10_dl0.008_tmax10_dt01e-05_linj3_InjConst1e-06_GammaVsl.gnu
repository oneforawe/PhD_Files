#filename: plot_2DKT_FPinject_vsl_lmax10_dl0.008_tmax10_dt01e-05_linj3_InjConst1e-06_GammaVsl.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_2DKT_FPinject_vsl_lmax10_dl0.008_tmax10_dt01e-05_linj3_InjConst1e-06_GammaVsl.ps'
#set size 1.0,1.1

set title "^4He Vortex-Pair Distribution {/Symbol-Oblique G} vs Pair Separation Parameter {/Times-Italic l}"

set ylabel "{/Symbol-Oblique G} (units...)"
set xlabel "{/Times-Italic l} (unitless)"
#set xrange [0:5]
#set yrange [1e-120:1e-40]
#set xrange [2:4]
#set yrange [1e-115:1e-75]
#set xrange [1.5:4.5]
#set yrange [1e-124:1e-62]
#set xrange [2.9:3.1]
#set yrange [1e-100:1e-90]
#set xrange [0:4]
#set yrange [1e-34:1]
##set xrange [2.75:3.25]
##set yrange [8e-8:1.05e-7]
set xrange [0:6]
########################
#set yrange [1e-17:1e-6]
#     set yrange [1e-40:1e-6]
set logscale y
set grid
set key bottom left
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot '2DKT_FPinject_vsl_IC0_T0.1_lmax5_dl0.004_tmax10000_dt01e-05_linj3_InjConst1e+200.dat' using 3:6 with linespoints title " T=0.10*Tc; InjConst=1.0e200 ", \
     '2DKT_FPinject_vsl_IC1_T0.1_lmax5_dl0.004_tmax10000_dt01e-05_linj3_InjConst1.7e-06_newBC.dat' using 3:6 with lines title " (newBC) T=0.10*Tc; InjConst=1.7e-06 ", \
     '2DKT_FPquench_vsl_T_1_0.1_lmax10_dl0.008_tmax10_dt01e-05_BC1_equil.dat' using 3:6 with linespoints title " T=1.00*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax5_dl0.004_linj3_InjConst100.dat' using 3:6 with linespoints title " T=0.10*Tc; InjConst=1.0e+02 "

#     '2DKT_FPinject_vsl_IC1_T0.1_lmax5_dl0.004_tmax10000_dt01e-05_linj3_InjConst1.7e-06.dat' using 3:6 with linespoints title " T=0.10*Tc; InjConst=1.7e-06 ", \
#
#'2DKT_FPinject_vsl_IC0_T0.1_lmax5_dl0.004_tmax10000_dt01e-05_linj3_InjConst1e+100.dat' using 3:6 with linespoints title " T=0.10*Tc "
#'2DKT_FPinject_vsl_IC0_T0.1_lmax10_dl0.008_tmax10000_dt01e-05_linj3_InjConst1e+100.dat' using 3:6 with linespoints title " T=0.10*Tc "
#
#     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax5_dl0.004_linj3_InjConst5.dat' using 3:6 with linespoints title " T=0.10*Tc; InjConst=5.0e+00 ", \
#     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax5_dl0.004_linj3_InjConst10.dat' using 3:6 with linespoints title " T=0.10*Tc; InjConst=1.0e+01 ", \
#     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax5_dl0.004_linj3_InjConst50.dat' using 3:6 with linespoints title " T=0.10*Tc; InjConst=5.0e+01 ", \
