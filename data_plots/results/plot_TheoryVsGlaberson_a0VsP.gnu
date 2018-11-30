#filename: plot_TheoryVsGlaberson_a0VsP.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_TheoryVsGlaberson_a0VsP.ps'

set title "Smallest Vortex-Loop Diameter {/Times-Italic a}_0({/Times-Italic P}) / {/Times-Italic a}_0(0) at Critical Temperature {/Times-Italic T}_{/Symbol-Oblique l}({/Times-Italic P}) \n vs Pressure {/Times-Italic P}"

set ylabel "{/Times-Italic a}_0({/Times-Italic P}) / {/Times-Italic a}_0(0) (unitless)"
set xlabel "{/Times-Italic P} (bar)"
set xrange [0:30]
#set yrange [0:0.008]
set grid
set key inside top left box width -4 height 0.5 spacing 1.4
set rmargin 2

plot 'Glaberson_a0VsP.dat' using ($1*(300./438)*(1./14.5037738007)):(1+$2*(0.3/350)):(1+$3*(0.3/350)):(1+$4*(0.3/350)) with yerrorbars pt 4 title " Glaberson's data ", \
     'a0VsP.dat'      using 1:6 with linespoints title " VLT calculations "
#     'a0VsP.dat' using 1:5 with points title " ... "
#     'a0vsPressTc.dat' using ($1*1.01325):3 title " T= " with linespoints