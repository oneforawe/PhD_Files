#filename: plot_ring3_C_KrK0versusTv.gnu
reset
set terminal png
set output "plot_ring3_C_KrK0versusTv.png"
set title "Renormalized coupling constant fraction (Kr/K0) \n vs Temperature-variable (1-T/Tc) \n with different length-scale cuts-off"
set ylabel "Kr/K0"
set xlabel "1-T/Tc"
set logscale xy
#set xrange [0.005:1]
set yrange [1e-05:10]
set grid
set key right bottom box title "Let tv = 1-T/Tc \n and A = 2.3885"
plot "ring3_Dexp100_C_beta_Cc0.40.dat" using 2:9 with linespoints title " Cc = 0.40 ", \
     "ring3_Dexp100_C_beta_Cc0.50.dat" using 2:9 with linespoints title " Cc = 0.50 ", \
     "ring3_Dexp100_C_beta_Cc0.60.dat" using 2:9 with linespoints title " Cc = 0.60 ", \
     "ring3_Dexp100_C_beta_Cc0.70.dat" using 2:9 with linespoints title " Cc = 0.70 ", \
     "ring3_Dexp100_C_beta_Cc0.80.dat" using 2:9 with linespoints title " Cc = 0.80 ", \
     "ring3_Dexp100_C_beta_Cc0.90.dat" using 2:9 with linespoints title " Cc = 0.90 ", \
     "ring3_Dexp100_C_beta_Cc1.00.dat" using 2:9 with linespoints title " Cc = 1.00 ", \
     "ring3_Dexp100_C_beta_Cc1.03.dat" using 2:9 with linespoints title " Cc = 1.03 ", \
     "KrK0versusT_powerlaw.dat" using 2:3 with linespoints title " A*tv^0.67168", \
     "ring1_Dexp100_A_Cc1.03.dat" using 2:9 with linespoints title " Cc = 1.03 ! ", \
     "ring1_Dexp100_B_Cc1.03.dat" using 2:9 with linespoints title " Cc = 1.03 ! ", \
     "ring3_Dexp100_A_beta_Cc1.03_reduced.dat" using 2:9 with linespoints title " Cc = 1.03 reduced"
