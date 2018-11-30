#filename: plot_2DKT_FPquench_adapt_lmax10_dl0.008_tmax10000_dt00.00018_GammaVsa.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_2DKT_FPquench_adapt_lmax10_dl0.008_tmax10000_dt00.00018_GammaVsa.ps'
#set size 1.0,1.1

set title "^4He Vortex-Pair Probability Density {/Symbol-Oblique G} vs Pair Separation {/Times-Italic a}/{/Times-Italic a}_0"

set ylabel "{/Symbol-Oblique G} (units...)"
set xlabel "{/Times-Italic a}/{/Times-Italic a}_0 (unitless)"
#set xrange [0:4e-5]
#set xrange [0:4.38e-4]
#set yrange [0:3e-5]
set logscale xy
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot '2DKT_FPquench_adapta_lmax10_dl0.008_tmax10000_dt00.00018.dat' using (exp($3)):6 with linespoints, \
     '2DKT_FPquench_adapt_lmax10_dl00.008_tmax10000_dt00.00018_largeDRE.dat' using (exp($3)):6 with linespoints, \
     '2DKT_FPquench_adapt_lmax10_dl00.008_tmax10000_dt00.00018_smallDRE.dat' using (exp($3)):6 with linespoints, \
     '2DKT_FPquench_adapt_lmax10_dl00.008_tmax10000_dt00.00018_DRE1e-15.dat' using (exp($3)):6 with linespoints, \
     '2DKT_FPquench_adapt_lmax10_dl02e-05_tmax10000_dt00.00018_DRE1e-15_TOL1e-26.dat' using (exp($3)):6 with linespoints, \
     '2DKT_FPquench_adapt_lmax10_dl02e-05_tmax10000_dt00.00018_DRE1e-15.dat' using (exp($3)):6 with linespoints, \
     '2DKT_FPquench_adapt_lmax10_dl02e-05_tmax10000_dt00.00018_DRE1e-15_TOL1e-05.dat' using (exp($5)):8 with linespoints

#'2DKT_FPquench_adapt_b4_lmax10_dl0.008_tmax10000_dt00.00018.dat' using (exp($3)):6 with linespoints, \

#     '2DKT_FPquench_adapt_lmax10_dl00.008_tmax10000_dt00.00018_DRE1e-10a.dat' using (exp($3)):6 with linespoints, \
#     '2DKT_FPquench_adapt_lmax10_dl00.008_tmax10000_dt00.00018_DRE1e-10.dat' using (exp($3)):6 with linespoints, \
#     '2DKT_FPquench_adapt_lmax10_dl00.008_tmax10000_dt00.00018_DRE1e-12.dat' using (exp($3)):6 with linespoints, \

#     '2DKT_FPquench_adaptb_lmax10_dl0.008_tmax10000_dt00.00018.dat' using (exp($3)):6 with linespoints, \
#     '2DKT_FPquench_adaptc_lmax10_dl0.008_tmax10000_dt00.00018.dat' using (exp($3)):6 with linespoints, \

# I think the adapt_b4 file is from the b4 program (2DKT_FPquench_adapt_b4.c)