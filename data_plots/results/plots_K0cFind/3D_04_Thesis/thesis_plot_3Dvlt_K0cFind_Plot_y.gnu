#filename: thesis_plot_3Dvlt_K0cFind_Plot_y.gnu
reset
set terminal postscript color eps enhanced
set output 'thesis_plot_3Dvlt_K0cFind_Plot_y.ps'
set size 0.85,0.65

set title "^4He Fugacity {/Times-Italic y} vs Length Scale {/Times-Italic l}"

set ylabel "{/Times-Italic y} (unitless)"
set xlabel "{/Times-Italic l} (unitless)"
#set xrange [0:4e-8]
set yrange [0:0.12]
set grid
set key inside top right box width -6 height 0.5 spacing 1.6 title " Finding {/Times-Italic K}@^{{} {/Times-Italic l}0^{^{}}}_{c  }, given a {/Times-Italic C}_c^ "
#set key at graph 0.86,-0.2 box width -5 height 0.5 spacing 1.2 title " Finding {/Times-Italic K}_{0c}, given a {/Times-Italic C}@_c^ "
#set bmargin 9

set pointsize 0.5

plot '3Dvlt_K0cFind_Plot_Cc1.10_K015.dat' using 1:3 with points title " {/Times-Italic C}_c = 1.1;   found {/Times-Italic K}@^{{} {/Times-Italic l}0^{^{}}}_{c } = 0.296468939007 ", \
     '3Dvlt_K0cFind_Plot_Cc0.50_K015.dat' using 1:3 with points title " {/Times-Italic C}_c = 0.5;   found {/Times-Italic K}@^{{} {/Times-Italic l}0^{^{}}}_{c } = 0.495108289977 "

# plot '3Dvlt_K0cFind_Plot_Cc1.10_K015.dat' using 1:3 with points title " {/Times-Italic C}_c = 1.1; found {/Times-Italic K}@^{{} {/Times-Italic l}0^{^{}}}_{c } = 0.296468939006784 ", \
#      '3Dvlt_K0cFind_Plot_Cc0.50_K015.dat' using 1:3 with points title " {/Times-Italic C}_c = 0.5; found {/Times-Italic K}@^{{} {/Times-Italic l}0^{^{}}}_{c } = 0.495108289977297 "
