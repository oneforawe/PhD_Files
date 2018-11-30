#filename: plot_RotonEnergyAtTcVsPressure.gnu
reset
set terminal png
set output 'plot_RotonEnergyAtTcVsPressure.png'
set title "Roton Energy Gap (Del/kB) at Tc vs Pressure"
set ylabel "Del/kB at Tc (K)"
set xlabel "P (atm)"
#set xrange [1:2.21]
set yrange [0:9]
set grid
set key left bottom box title " Del/kB at various values of Tc "
plot 'DelokVsPressTc.dat' using 1:3 title " Calculated data from RhosAmps.ods " with linespoints, \
     'RotonEnergyVsPressTc.dat' using 1:3 title " Data taken/extrapolated from Brooks, Donnelly" with linespoints