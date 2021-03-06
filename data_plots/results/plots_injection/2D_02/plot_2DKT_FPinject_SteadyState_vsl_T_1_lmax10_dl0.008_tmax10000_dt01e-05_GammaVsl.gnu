#filename: plot_2DKT_FPinject_SteadyState_vsl_T_1_lmax10_dl0.008_tmax10000_dt01e-05_GammaVsl.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_2DKT_FPinject_SteadyState_vsl_T_1_lmax10_dl0.008_tmax10000_dt01e-05_GammaVsl.ps'
#set size 1.0,1.1

set title "^4He Vortex-Pair Distribution {/Symbol-Oblique G} vs Pair Separation Parameter {/Times-Italic l}"

set ylabel "{/Symbol-Oblique G} (units...)"
set xlabel "{/Times-Italic l} (unitless)"
#set xrange [0:5]
#set yrange [1e-8:10]
set yrange [1e-40:1e4]
#set yrange [1e-100:1e4]
#set yrange [1e-40:10]
set logscale y
set grid
set key outside center right box width -5
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot '2DKT_FPquench_vsl_T_0.1_0.1_lmax10_dl0.008_tmax10_dt01e-05_BC1_equil.dat' using 3:6 with lines title " no injection; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst1e-35.dat' using 3:6 with lines title " InjConst=1.0e-35; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst1e-33.dat' using 3:6 with lines title " InjConst=1.0e-33; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst1e-26.dat' using 3:6 with lines title " InjConst=1.0e-26; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst1e-20.dat' using 3:6 with lines title " InjConst=1.0e-20; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst1e-13.dat' using 3:6 with lines title " InjConst=1.0e-13; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst1e-07.dat' using 3:6 with lines title " InjConst=1.0e-07; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst1e-06.dat' using 3:6 with lines title " InjConst=1.0e-06; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst1.5e-06.dat' using 3:6 with lines title " InjConst=1.5e-06; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst1.7e-06.dat' using 3:6 with lines title " InjConst=1.7e-06; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst1.8e-06.dat' using 3:6 with lines title " InjConst=1.8e-06; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst1.9e-06.dat' using 3:6 with lines title " InjConst=1.9e-06; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst2e-06.dat' using 3:6 with lines title " InjConst=2.0e-06; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst2.1e-06.dat' using 3:6 with lines title " InjConst=2.1e-06; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst2.5e-06.dat' using 3:6 with lines title " InjConst=2.5e-06; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst3e-06.dat' using 3:6 with lines title " InjConst=3.0e-06; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst5e-06.dat' using 3:6 with lines title " InjConst=5.0e-06; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst1e-05.dat' using 3:6 with lines title " InjConst=1.0e-05; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst5e-05.dat' using 3:6 with lines title " InjConst=5.0e-05; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst0.0001.dat' using 3:6 with lines title "InjConst=1.0e-04; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst0.001.dat' using 3:6 with lines title " InjConst=1.0e-03; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst0.01.dat' using 3:6 with lines title " InjConst=1.0e-02; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst0.05.dat' using 3:6 with lines title " InjConst=5.0e-02; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst0.1.dat' using 3:6 with lines title " InjConst=1.0e-01; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst0.5.dat' using 3:6 with lines title " InjConst=5.0e-01; T=0.10*Tc ", \
     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst1.dat' using 3:6 with lines title " InjConst=1; T=0.10*Tc ", \
     '2DKT_FPquench_vsl_T_1_0.1_lmax10_dl0.008_tmax10_dt01e-05_BC1_equil.dat' using 3:6 with lines title " no injection; T=1.00*Tc "

#     '2DKT_FPinject_SteadyState_vsl_T_1_lmax10_dl0.008_linj3_InjConst1e-06.dat' using 3:6 with lines title " T=1.00*Tc (InjConst=1e-6) ", \
#     '2DKT_FPinject_SteadyState_vsl_T_1_lmax10_dl0.008_linj3_InjConst1.dat' using 3:6 with lines title " T=1.00*Tc (InjConst=1) "
#
#     '2DKT_FPinject_SteadyState_vsl_T_1_lmax10_dl0.008_linj3_InjConst1e-06_sinkzf.dat' using 3:6 with lines title " T=1.00*Tc (InjConst=1e-6) w/ sinkzf"
#     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst1e-20_sinkzf.dat' using 3:6 with lines title " T=0.10*Tc (InjConst=1e-20) w/ sinkzf", \
#
#     '2DKT_FPquench_inject_SteadyState_vsl_T_1_lmax10_dl0.008_linj3_InjConst1.dat' using 3:6 with lines title " T=1.00*Tc (InjConst=1) ", \
#     '2DKT_FPquench_inject_SteadyState_vsl_T_1_lmax10_dl0.008_linj3_InjConst1e-06.dat' using 3:6 with lines title " T=1.00*Tc (InjConst=1e-6) ", \
#     '2DKT_FPinject_SteadyState_vsl_T_1_lmax10_dl0.008_linj3_InjConst0.0001.dat' using 3:6 with lines title " T=1.00*Tc (InjConst=1e-4) ", \
#     '2DKT_FPinject_SteadyState_vsl_T_1_lmax10_dl0.008_linj3_InjConst0.001.dat' using 3:6 with lines title " T=1.00*Tc (InjConst=1e-3) ", \
#     '2DKT_FPinject_SteadyState_vsl_T_1_lmax10_dl0.008_linj3_InjConst1e-06_sink.dat' using 3:6 with lines title " T=1.00*Tc (InjConst=1e-4) w/ sink"
#     '2DKT_FPinject_SteadyState_vsl_T_1_lmax10_dl0.008_linj3_InjConst1e-06_sink.dat' using 3:6 with lines title " T=1.00*Tc (InjConst=1e-6) w/ sink", \
#     '2DKT_FPinject_SteadyState_vsl_T_0.1_lmax10_dl0.008_linj3_InjConst1e-20_sink.dat' using 3:6 with lines title " T=0.10*Tc (InjConst=1e-20) w/ sink", \
