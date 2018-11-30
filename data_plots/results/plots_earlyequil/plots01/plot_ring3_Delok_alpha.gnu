#filename: plot_ring3_Delok_alpha.gnu
reset
set terminal png
set output 'plot_ring3_Delok_alpha.png'
set title "Roton Energy Gap (Del/kB) \n vs Temperature (T/Tc)"
set ylabel "Del/kB (K)"
set xlabel "T/Tc"
#set xrange [1:2.21]
set yrange [6.5:11.5]
set grid
set key left bottom box title "Fit of Donnelly data multiplied by ~6.8/5.3"
plot 'ring3_Dexp100_C_alpha.dat' using 1:11 title " computer fit " with linespoints
#, \
#     'funkyfactorfit.dat' using 1:2 title " factor" with linespoints
