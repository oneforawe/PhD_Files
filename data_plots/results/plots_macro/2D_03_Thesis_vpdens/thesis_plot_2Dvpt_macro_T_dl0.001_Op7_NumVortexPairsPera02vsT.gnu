#filename: thesis_plot_2Dvpt_macro_T_dl0.001_Op7_NumVortexPairsPera02vsT.gnu

# This is meant to show how the number of pairs increases with increasing temperature.

reset
set terminal postscript color eps enhanced
set output 'thesis_plot_2Dvpt_macro_T_dl0.001_Op7_NumVortexPairsPera02vsT.ps'
set size 0.85,0.85

set title "^4He Number of Vortex Pairs per {/Times-Italic a}_0^2 vs Temperature Fraction {/Times-Italic T} / {/Times-Italic T}_{KT}"

set ylabel "{/Times-Italic (N / A)} {{/Times-Italic a}_0}^2 (unitless)"
set xlabel "{/Times-Italic T} / {/Times-Italic T}_{KT} (unitless)"
#set xrange [0:1]
#set yrange [-0.5:12]
##set yrange [-0.1:2]
set yrange [-0.5e-1*0.1:0.5e-1*2]
##set yrange [-0.05*1:1]
#set logscale y
#set grid
#set key inside left bottom box width -2 height 0.5 spacing 1.2 title " ^{}Vortex Pair Theory "
set key at graph 0.34,0.48 box width -2 height 0.5 spacing 1.2 title " ^{}Vortex Pair Theory "
#set rmargin 35

#set pointsize 0.5

plot '2Dvpt_macro_T_lmax0_dl-nan_Op7.dat'    using 1:(-$8) with lines title " {/Times-Italic l}_{max} = 0 ", \
     '2Dvpt_macro_T_lmax1_dl0.001_Op7.dat'   using 1:(-$8) with lines title " {/Times-Italic l}_{max} = 1 ", \
     '2Dvpt_macro_T_lmax2_dl0.001_Op7.dat'   using 1:(-$8) with lines title " {/Times-Italic l}_{max} = 2 ", \
     '2Dvpt_macro_T_lmax4_dl0.001_Op7.dat'   using 1:(-$8) with lines title " {/Times-Italic l}_{max} = 4 ", \
     '2Dvpt_macro_T_lmax10_dl0.001_Op7.dat'  using 1:(-$8) with lines title " {/Times-Italic l}_{max} = 10 ", \
     '2Dvpt_macro_T_lmax83_dl0.001_Op7.dat'  using 1:(-$8) with lines title " {/Times-Italic l}_{max} = 83 "

#     '3Dvlt_macro_T_Cc1.11_lmax100_dl0.001_Op7.dat' using 1:(-(4.0/pi)**2*$8) with lines title " {/Times-Italic l}_{max} = 100 "
