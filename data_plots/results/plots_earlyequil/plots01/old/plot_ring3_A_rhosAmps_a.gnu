#filename: plot_ring3_A_rhosAmps_a.gnu
reset
set terminal png
set output "plot_ring3_A_rhosAmps_a.png"
set title "Rho_s Amplitudes (A' = {rhos/rho}/tv^nu = {Kr/K0}/tv^nu) \n vs temperature variable (tv = 1-T/Tc) \n with various critical lambda core energies (Cc)"
set ylabel "rhos amp"
set xlabel "tv"
#set logscale xy
#set xrange [0:0.01]
#set yrange [0:30]
set grid
set key right outside box title "nu = 0.6716883"
plot "ring3_Dexp100_A_beta_Cc1.03.dat" using 2:($9/($2**0.6716883)) with linespoints title " Cc = 1.03 ", \
     "ring3_Dexp100_A_beta_Cc1.00.dat" using 2:($9/($2**0.6716883)) with linespoints title " Cc = 1.00 ", \
     "ring3_Dexp100_A_beta_Cc0.90.dat" using 2:($9/($2**0.6716883)) with linespoints title " Cc = 0.90 ", \
     "ring3_Dexp100_A_beta_Cc0.80.dat" using 2:($9/($2**0.6716883)) with linespoints title " Cc = 0.80 ", \
     "ring3_Dexp100_A_beta_Cc0.70.dat" using 2:($9/($2**0.6716883)) with linespoints title " Cc = 0.70 ", \
     "ring3_Dexp100_A_beta_Cc0.60.dat" using 2:($9/($2**0.6716883)) with linespoints title " Cc = 0.60 ", \
     "ring3_Dexp100_A_beta_Cc0.50.dat" using 2:($9/($2**0.6716883)) with linespoints title " Cc = 0.50 ", \
     "ring3_Dexp100_A_beta_Cc0.40.dat" using 2:($9/($2**0.6716883)) with linespoints title " Cc = 0.40 "
     