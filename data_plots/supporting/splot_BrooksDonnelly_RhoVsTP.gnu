#filename: splot_BrooksDonnelly_RhoVsTP.gnu
reset
set terminal postscript color eps enhanced
set output 'splot_BrooksDonnelly_RhoVsTP.ps'
#set size 1.3,1

set title "^4He Density {/Symbol-Oblique r} \n vs Temperature {/Times-Italic T} and Pressure {/Times-Italic P} "

set zlabel "{/Symbol-Oblique r} (kg m^{-3})"
set xlabel "{/Times-Italic P} (bar)"
set ylabel "{/Times-Italic T} (K)"
#set logscale x
set zrange [130:180]
#set grid
set key top box spacing 1.2 title " Data from Brooks/Donnelly "
set lmargin 8

set pm3d
splot 'BrooksDonnelly_RhoVsTP_3d.dat' with pm3d
