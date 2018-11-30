#filename: plot_2Dvpt_quench_vsl_T_1_0.1_lmax10_dl0.001_tmax10_dt2.5e-07_BC1_GammaVsl.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_2Dvpt_quench_vsl_T_1_0.1_lmax10_dl0.001_tmax10_dt2.5e-07_BC1_GammaVsl.ps'
#set size 1.0,1.1

set title "^4He Vortex-Pair Distribution {/Symbol-Oblique G} vs Pair Separation Parameter {/Times-Italic l}"

set ylabel "{/Symbol-Oblique G} (units...)"
set xlabel "{/Times-Italic l} (unitless)"
#set xrange [0:2]
#set yrange [1e-10:1e-3]
set logscale y
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot '2Dvpt_quench_vsl_T_1_0.1_lmax10_dl0.001_tmax10_dt2.5e-07_BC1.dat' using 3:6 with linespoints title " {/Symbol-Oblique G}, BC1 1.00 to 0.1 of Tc ", \
     '2Dvpt_quench_vsl_T_1_0.1_lmax10_dl0.001_tmax10_dt2.5e-07_BC1_t0.dat' using 3:(exp(-6.283185*$4*$3)) with linespoints title " exp(-2*pi*{/Times-Italic K}_i*{/Times-Italic l}) ", \
     '2Dvpt_quench_vsl_T_1_0.1_lmax10_dl0.001_tmax10_dt2.5e-07_BC1_t0.dat' using 3:($6/(exp(-6.283185*$4*$3))) with linespoints title " {/Symbol-Oblique G}_i/exp(-2*pi*{/Times-Italic K}_i*{/Times-Italic l}) ", \
     exp(-6*x) with lines title " exp(-6*{/Times-Italic l}) "