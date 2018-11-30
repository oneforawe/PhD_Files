#filename: plot_3Dvlt_quench_vsl_T_1_0.7_P_00.050_lmax10_dl0.002_tmax100_dt2e-06_BC1_GammaVsl.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_3Dvlt_quench_vsl_T_1_0.7_P_00.050_lmax10_dl0.002_tmax100_dt2e-06_BC1_GammaVsl.ps'
#set size 1.0,1.1

set title "^4He Vortex-Loop Distribution {/Symbol-Oblique G} vs Loop Diameter Scale {/Times-Italic l}"

set ylabel "{/Symbol-Oblique G} (vortex-loops/{/Times-Italic a}_0^6)"
set xlabel "{/Times-Italic l} (unitless)"
#set xrange [0:1]
#set yrange [1e-30:1]
set logscale y
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot '3Dvlt_quench_vsl_T_1_0.7_P_00.050_lmax10_dl0.002_tmax100_dt2e-06_BC1.dat'  using 3:7 with linespoints  title " BC1 1.00 to 0.7 of Tc ", \
     '3Dvlt_quench_vsl_T_0.7_0.7_P_00.050_lmax10_dl0.002_tmax100_dt2e-06_BC1.dat'  using 3:7 with linespoints  title " new equil at 0.7 of Tc "
