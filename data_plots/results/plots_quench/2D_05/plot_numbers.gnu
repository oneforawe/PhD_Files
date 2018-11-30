#filename: plot_numbers.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_numbers.ps'
#set size 1.0,1.1

set title "^4He Vortex-Pair Distribution {/Symbol-Oblique G} vs Pair Separation Parameter {/Times-Italic l}"

set ylabel "{/Symbol-Oblique G} (units...)"
set xlabel "{/Times-Italic l} (unitless)"
set xrange [0:2]
#set yrange [1e-10:1e-3]
#set xrange [0:0.5]
#set yrange [1e-4:1e-3]
# wider look
#set xrange [0:1]
#set yrange [1e-6:1e-3]
# yet wider
##set xrange [0:2]
##set yrange [1e-8:1e-3]
# different look
#set xrange [2:2.5]
#set yrange [1e-9:1e-7]
# different look
#set xrange [0.5:1.5]
#set yrange [1e-7:1e-4]
set logscale y
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot '2DKT_FPquench_vsl_T_1_0.5_lmax10_dl0.008_tmax10_dt01e-05_BC1.dat' using 3:6 with linespoints title " BC1 1.00 to 0.5 of Tc ", \
     '2DKT_FPquench_vsl_T_1_0.5_lmax10_dl0.008_tmax10_dt01e-05_BC2.dat' using 3:6 with linespoints title " BC2 1.00 to 0.5 of Tc ", \
     '2DKT_FPquench_vsl_T_1_0.5_lmax10_dl0.008_tmax10_dt01e-05_BC3.dat' using 3:6 with linespoints title " BC3 1.00 to 0.5 of Tc ", \
     1.0e-06 with lines, \
     2.0e-06 with lines, \
     5.0e-06 with lines, \
     8.0e-06 with lines, \
     9.0e-06 with lines, \
     9.9e-06 with lines, \
     1.0e-05 with lines

#     '2DKT_FPquench_vsl_T_1_0.5_lmax10_dl0.008_tmax10_dt01e-05_BC2_other.dat' using 3:6 with linespoints title " BC2 1.00 to 0.5 of Tc ", \
#     '2DKT_FPquench_vsl_T_1_0.5_lmax10_dl0.008_tmax10_dt01e-05_BC3_wrong.dat' using 3:6 with linespoints title " BC3 1.00 to 0.5 of Tc "
#     '2DKT_FPquench_vsl_T_1_0.5_lmax10_dl0.008_tmax10_dt01e-05_BC3_orig1.dat' using 3:6 with linespoints title " BC3 1.00 to 0.5 of Tc ", \
#     '2DKT_FPquench_vsl_T_1_0.5_lmax10_dl0.008_tmax10_dt01e-05_BC3_orig2.dat' using 3:6 with linespoints title " BC3 1.00 to 0.5 of Tc ", \
#     '2DKT_FPquench_vsl_T_1_0.5_lmax10_dl0.008_tmax10_dt01e-05_BC3_orig3.dat' using 3:6 with linespoints title " BC3 1.00 to 0.5 of Tc ", \
# BC3_orig1: contradiction* in formulas, prob-dens allowed to fall below new equilibrium value
# BC3_orig2: contradiction* in formulas, prob-dens stopped at new equilibrium value
# BC3_orig3: no contradiction* (but perhaps an overestimate in prob-dens due to += G[0]*exp(2*l[0])*dl; )
# BC3_orig4: 0.5 was not in += 0.5*G[0]*exp(2*l[0])*dl; for equilibrium calculation (as with orig3,2,1)
# * 0.5 was in += 0.5*G[0]*exp(2*l[0])*dl; but not in (2*PI*K[0]*dl*dt + 0.5*t0*dl*dl - dt)
