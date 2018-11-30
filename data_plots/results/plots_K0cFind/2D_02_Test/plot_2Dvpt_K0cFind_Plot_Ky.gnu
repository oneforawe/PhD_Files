#filename: plot_2Dvpt_K0cFind_Plot_Ky.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_2Dvpt_K0cFind_Plot_Ky.ps'
set size 0.85,0.85

set title "^4He Film Coupling Parameter {/Times-Italic K} vs Fugacity {/Times-Italic y}"

set ylabel "{/Times-Italic K} (unitless)"
set xlabel "{/Times-Italic y} (unitless)"
set xrange [2.0/pi-1e-1:2.0/pi+1e-1]
set yrange [0:2e-4]
set grid
set key at graph 0.775,-0.22 box width -5 height 0.5 spacing 1.2 title " Finding {/Times-Italic K}_{0c}, using the Villain model^{} "
#set key outside below center box width -5 height 0.5 spacing 1.2 title " Finding {/Times-Italic K}_{0c}, using the Villain model^{} "
set bmargin 9

set pointsize 0.5

plot '2Dvpt_K0cFind_Plot_Output.dat' using 2:3 with points title " Found {/Times-Italic K}_{0c} = 0.7478523792 "
