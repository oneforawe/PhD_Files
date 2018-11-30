#filename: splot_a0VsTP.gnu
reset
set terminal postscript color eps enhanced
set output 'splot_a0VsTP.ps'
#set size 1.3,1

set title "^4He Minimal Vortex Diameter {/Times-Italic a}_0 vs Temperature {/Times-Italic T} and Pressure {/Times-Italic P} "

set zlabel "{/Times-Italic a}_0 (angstrom)"
set xlabel "{/Times-Italic P} (bar)"
set ylabel "{/Times-Italic T} (K)"
#set logscale x
#set xrange [0:5]
#set grid
set key outside right box width -10 spacing 1.2 title " Data from {/CM-Typewriter RhosAmps.ods} "
#set rmargin 17

set pm3d
splot 'a0VsTP_3d.dat' with pm3d
