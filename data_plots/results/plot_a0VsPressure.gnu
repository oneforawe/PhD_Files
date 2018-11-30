#filename: plot_a0VsPressure.gnu
reset
set terminal png
set output 'plot_a0VsPressure.png'
set title "Smallest Vortex-Ring Diameter (a0) vs Pressure"
set ylabel "a0 (?)"
set xlabel "P (atm)"
#set xrange [1:2.21]
set yrange [0:9]
set grid
set key left bottom box title "Data from Brooks, Donnelly"
plot 'a0VsPressTc.dat' using 1:3 title " T= " with linespoints, \
     'a0VsPressTc.dat' using 1:3 title " T= " with linespoints, \
