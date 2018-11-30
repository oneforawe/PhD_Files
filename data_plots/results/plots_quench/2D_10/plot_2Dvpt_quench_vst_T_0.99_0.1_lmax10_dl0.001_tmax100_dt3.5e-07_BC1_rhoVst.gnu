#filename: plot_2Dvpt_quench_vst_T_0.99_0.1_lmax10_dl0.001_tmax100_dt3.5e-07_BC1_rhoVst.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_2Dvpt_quench_vst_T_0.99_0.1_lmax10_dl0.001_tmax100_dt3.5e-07_BC1_rhoVst.ps'
#set size 1.0,1.1

set title "^4He Vortex-Pair Density {/Symbol-Oblique r} vs Time {/Symbol-Oblique t}"

set ylabel "{/Symbol-Oblique r} (units...)"
set xlabel "{/Symbol-Oblique t} (unitless)"
set xrange [1e-7:10000]
#set yrange [1e-10:1e-3]
set logscale xy
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot '2Dvpt_quench_vst_T_0.99_0.1_lmax10_dl0.001_tmax100_dt3.5e-07_BC1.dat' using 2:9 with linespoints title " {/Symbol-Oblique r}, BC1 0.99 to 0.1 of Tc ", \
     '2Dvpt_quench_vst_T_0.99_0.1_lmax10_dl0.001_tmax100_dt3.5e-07_BC1.dat' using 2:(0.00001*$2**$12) with linespoints title " 10^{-5}*({/Times-Italic t}/{/Symbol-Oblique t})^{(-zscale/z)} ", \
     0.00001/(x*log(x)) with lines title " 10^{-5}/({/Symbol-Oblique t}*ln({/Symbol-Oblique t})) ", \
     0.0001/x with lines title " 10^{-4}/{/Symbol-Oblique t} ", \
     0.0001/(log(x)) with lines title " 10^{-4}/(ln({/Symbol-Oblique t})) "
