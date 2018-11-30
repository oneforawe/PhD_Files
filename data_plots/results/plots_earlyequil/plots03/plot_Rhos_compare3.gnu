#filename: plot_Rhos_compare3.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_Rhos_compare3.ps'

set title "Superfluid fraction {/Symbol-Oblique r}_s/{/Symbol-Oblique r}_s^0 \n or Renormalized coupling constant fraction ({/Times-Italic K}_r/{/Times-Italic K}_0) \n vs Temperature {/Times-Italic T}"

set ylabel "{/Symbol-Oblique r}_s/{/Symbol-Oblique r}_s^0 = {/Times-Italic K}_r/{/Times-Italic K}_0 (unitless)"
set xlabel "{/Times-Italic T} (K)"
set grid
set key left bottom box

plot 'Donnelly_rhos.dat' using 1:2 with points title " Donnelly data  ", \
     'vlt_HeatCap_P_00.050_DK0_5e-07_i47.dat' using (2.1768*$1):8 with linespoints title " new calculation (un/corrected?) "
#     'vlt_HeatCap_P_00.050_DK0_5e-07_i47.dat' using (2.172*$1):8 with linespoints title " new calculation (un/corrected?) "
#     'vlt_HeatCap_P_00.050_DK0_5e-07_i47.dat' using (2.1768*$1):8 with linespoints title " new calculation (un/corrected?) "
