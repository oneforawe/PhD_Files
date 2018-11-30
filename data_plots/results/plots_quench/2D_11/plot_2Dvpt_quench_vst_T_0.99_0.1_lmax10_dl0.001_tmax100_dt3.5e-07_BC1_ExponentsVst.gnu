#filename: plot_2Dvpt_quench_vst_T_0.99_0.1_lmax10_dl0.001_tmax100_dt3.5e-07_BC1_ExponentsVst.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_2Dvpt_quench_vst_T_0.99_0.1_lmax10_dl0.001_tmax100_dt3.5e-07_BC1_ExponentsVst.ps'
#set size 1.0,1.1

set title "^4He Vortex-Pair Density Exponents (theory and calculated) vs Time {/Times-Italic t}/{/Symbol-Oblique t}"

set ylabel "Exponents (units...)"
set xlabel "{/Times-Italic t}/{/Symbol-Oblique t} (unitless)"
#set xrange [0:2]
#set yrange [1e-10:1e-3]
set logscale x
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot '2Dvpt_quench_vst_T_0.99_0.1_lmax10_dl0.001_tmax100_dt3.5e-07_BC1_new.dat' using 2:10 with linespoints title " calculated = d(ln {/Symbol-Oblique r})/d(ln {/Times-Italic t}) ", \
     '2Dvpt_quench_vst_T_0.99_0.1_lmax10_dl0.001_tmax100_dt3.5e-07_BC1_new.dat' using 2:13 with linespoints title " 1a theory = -zscale/z ", \
     '2Dvpt_quench_vst_T_0.99_0.1_lmax10_dl0.001_tmax100_dt3.5e-07_BC1_new.dat' using 2:14 with linespoints title " 1b theory = -zscale/z ", \
     '2Dvpt_quench_vst_T_0.99_0.1_lmax10_dl0.001_tmax100_dt3.5e-07_BC1_new.dat' using 2:15 with linespoints title " 2a theory = -zscale/z ", \
     '2Dvpt_quench_vst_T_0.99_0.1_lmax10_dl0.001_tmax100_dt3.5e-07_BC1_new.dat' using 2:16 with linespoints title " 2b theory = -zscale/z "
