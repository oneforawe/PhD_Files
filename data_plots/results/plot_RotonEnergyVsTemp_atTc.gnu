#filename: plot_RotonEnergyVsTemp_atTc.gnu
reset
set terminal png
set output 'plot_RotonEnergyVsTemp_atTc.png'
set title "Roton Energy Gap (Del/kB) vs Temperature (Tc)"
set ylabel "Del/kB (K)"
set xlabel "T (K)"
#set xrange [1:2.21]
set yrange [0:9]
set grid
set key left bottom box title "Curves from Brooks, Donnelly \n Delok points from RhosAmps.ods"
plot 'RotonEnergyVsTemp_00atm.txt' using 1:2 title " 00 atm " with linespoints, \
     'RotonEnergyVsTemp_01atm.txt' using 1:2 title " 01 atm " with linespoints, \
     'RotonEnergyVsTemp_05atm.txt' using 1:2 title " 05 atm " with linespoints, \
     'RotonEnergyVsTemp_10atm.txt' using 1:2 title " 10 atm " with linespoints, \
     'RotonEnergyVsTemp_15atm.txt' using 1:2 title " 15 atm " with linespoints, \
     'RotonEnergyVsTemp_20atm.txt' using 1:2 title " 20 atm " with linespoints, \
     'RotonEnergyVsTemp_24atm.txt' using 1:2 title " 24 atm " with linespoints, \
     'RotonEnergyVsPressTc.dat' using 2:3 title " values at Tc " with points 10 5, \
     'DelokVsPressTc.dat' using 2:3 title " Delok at Tc " with points
#     'DelokVsTc.dat' using 2:($3-1.5) title " Delok at Tc - 1.5" with points
