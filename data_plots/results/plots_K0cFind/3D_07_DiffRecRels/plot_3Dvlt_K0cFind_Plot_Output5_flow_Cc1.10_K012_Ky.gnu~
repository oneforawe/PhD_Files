#filename: thesis_plot_3Dvlt_K0cFind_Plot_Ky.gnu
reset
set terminal postscript color eps enhanced
set output 'thesis_plot_3Dvlt_K0cFind_Plot_Ky.ps'
set size 0.85,0.85

set title "^4He Fugacity {/Times-Italic y} vs Coupling Parameter {/Times-Italic K}"

set ylabel "{/Times-Italic y} (unitless)"
set xlabel "{/Times-Italic K} (unitless)"
set xrange [0.38750818971-3e-4:0.38750818971+3e-4]
set yrange [0.062421005458-5e-5:0.062421005458+5e-5] 
set grid
set key at graph 0.85,-0.21 box width -5 height 0.5 spacing 1.6 title " Finding {/Times-Italic K}@^{{} {/Times-Italic l}0^{^{}}}_{c  }, given a {/Times-Italic C}_c^ "
##set key at graph 0.875,-0.21 box width -5 height 0.5 spacing 1.6 title " Finding {/Times-Italic K}@^{{} {/Times-Italic l}0^{^{}}}_{c  }, given a {/Times-Italic C}_c^ "
#set key outside below center box width -5 height 0.5 spacing 1.2 title " Finding {/Times-Italic K}_{0c}, given a {/Times-Italic C}@_c^ "
set bmargin 9.7

set pointsize 0.5

plot '3Dvlt_K0cFind_Plot_Cc1.10_K015.dat' using 2:3 with points title " {/Times-Italic C}_c = 1.1;   found {/Times-Italic K}@^{{} {/Times-Italic l}0^{^{}}}_{c } = 0.296468939007 ", \
     '3Dvlt_K0cFind_Plot_Cc0.50_K015.dat' using 2:3 with points title " {/Times-Italic C}_c = 0.5;   found {/Times-Italic K}@^{{} {/Times-Italic l}0^{^{}}}_{c } = 0.495108289977 "

# plot '3Dvlt_K0cFind_Plot_Cc1.10_K015.dat' using 2:3 with points title " {/Times-Italic C}_c = 1.1; we find {/Times-Italic K}@^{{} {/Times-Italic l}0^{^{}}}_{c } = 0.296468939006784 ", \
#      '3Dvlt_K0cFind_Plot_Cc0.50_K015.dat' using 2:3 with points title " {/Times-Italic C}_c = 0.5; we find {/Times-Italic K}@^{{} {/Times-Italic l}0^{^{}}}_{c } = 0.495108289977297 "
