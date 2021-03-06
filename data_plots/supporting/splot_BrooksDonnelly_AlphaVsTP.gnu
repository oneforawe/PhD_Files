#filename: splot_BrooksDonnelly_AlphaVsTP.gnu
reset
set terminal postscript color eps enhanced
set output 'splot_BrooksDonnelly_AlphaVsTP.ps'
#set size 1.3,1

set title "^4He Thermal Expansion Coefficient {/Symbol-Oblique a} \n vs Temperature {/Times-Italic T} and Pressure {/Times-Italic P} "

set zlabel "{/Symbol-Oblique a} (K^{-1})"
set xlabel "{/Times-Italic P} (bar)"
set ylabel "{/Times-Italic T} (K)"
#set logscale x
set zrange [-0.1:0.01]
#set grid
set key top box spacing 1.2 title " Data from Brooks/Donnelly "
set lmargin 5

set pm3d
splot 'BrooksDonnelly_AlphaVsTP_3d.dat' with pm3d
