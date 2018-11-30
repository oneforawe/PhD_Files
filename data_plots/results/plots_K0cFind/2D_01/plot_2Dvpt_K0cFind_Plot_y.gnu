#filename: plot_2Dvpt_K0cFind_Plot_y.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_2Dvpt_K0cFind_Plot_y.ps'
set size 0.85,0.65

set title "^4He Film Fugacity {/Times-Italic y} vs Length Scale {/Times-Italic l}"

set ylabel "{/Times-Italic y} (unitless)"
set xlabel "{/Times-Italic l} (unitless)"
#set xrange [0:4e-8]
set yrange [1e-15:1e5]
set logscale y
set grid
set key inside top right box width -6 height 0.5 spacing 1.2 title " Finding {/Times-Italic K}_{0c}, using the Villain model^{} "
#set key at graph 0.86,-0.2 box width -5 height 0.5 spacing 1.2 title " Finding {/Times-Italic K}_{0c}, using the Villain model^{} "
#set bmargin 9

set pointsize 0.5

plot '2Dvpt_K0cFind_Plot_Output.dat' using 1:3 with points title " Found {/Times-Italic K}_{0c} = 0.74785238 "
