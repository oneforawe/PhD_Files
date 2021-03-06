#filename: plot_3D_KrK0versusTv_B_findA.gnu
reset
set terminal png
set output "plot_3D_KrK0versusTv_B_findA.png"
set title "Renormalized coupling constant fraction (Kr/K0) \n divided by (1-T/Tc)^0.6716883 \n vs Temperature-variable (1-T/Tc) in 3D"
set ylabel "Kr/K0"
set xlabel "1-T/Tc"
set logscale x
set grid
set key left bottom box
plot "ring1_B_Dexp100.dat" using 2:($9/(($2)**(0.6716883))) with linespoints title " l_max = 100 "
