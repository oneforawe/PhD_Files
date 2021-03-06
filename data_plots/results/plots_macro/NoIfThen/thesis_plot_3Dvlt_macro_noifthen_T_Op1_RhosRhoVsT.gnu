#filename: thesis_plot_3Dvlt_macro_T_Cc1.11_dl0.001_Op7_RhosRhoVsT.gnu
reset
set terminal postscript color eps enhanced
set output 'thesis_plot_3Dvlt_macro_T_Cc1.11_dl0.001_Op7_RhosRhoVsT.ps'
set size 0.85,0.85

set title "^4He Superfluid Fraction {/Symbol-Oblique r}_s / {/Symbol-Oblique r} vs Temperature Fraction {/Times-Italic T} / {/Times-Italic T}_{/Symbol-Oblique l}"

set ylabel "{/Symbol-Oblique r}_s / {/Symbol-Oblique r} (unitless)"
set xlabel "{/Times-Italic T} / {/Times-Italic T}_{/Symbol-Oblique l} (unitless)"
#set xrange [0:1]
##set yrange [0:1.1]
set yrange [-0.1:1.1]
#set logscale y
#set grid
#set key inside left bottom box width -2 height 0.5 spacing 1.2 title " ^{}Vortex Loop Theory "
set key at graph 0.32,0.32 box width -2 height 0.5 spacing 1.2 title " ^{}Vortex Loop Theory "
#set rmargin 35

#set pointsize 0.5

plot '3Dvlt_macro_nitKth_T_Cc0.40_lmax10_dl0.001_Op1_usingK0cFind1.dat'     using 1:8 with lines lt 1 lw 3 lc rgb "yellow" title " {/Times-Italic l}_{max} = 10 ", \
     '3Dvlt_macro_nitKth_T_Cc0.40_lmax10_dl0.001_Op1_usingK0cFind3dl.dat'   using 1:8 with lines lt 2 lw 2 lc rgb "green"  title " {/Times-Italic l}_{max} = 10 ", \
     '3Dvlt_macro_nitKth_T_Cc0.40_lmax10_dl0.001_Op1_usingK0cFind3ds.dat'   using 1:8 with lines lt 3 lw 1 lc rgb "red"    title " {/Times-Italic l}_{max} = 10 ", \
     '3Dvlt_macro_nitKth_T_Cc1.10_lmax10_dl0.001_Op1_usingK0cFind1.dat'     using 1:8 with lines lt 1 lw 3 lc rgb "cyan"   title " {/Times-Italic l}_{max} = 10 ", \
     '3Dvlt_macro_nitKth_T_Cc1.10_lmax10_dl0.001_Op1_usingK0cFind3dl.dat'   using 1:8 with lines lt 2 lw 2 lc rgb "orange" title " {/Times-Italic l}_{max} = 10 ", \
     '3Dvlt_macro_nitKth_T_Cc1.10_lmax10_dl0.001_Op1_usingK0cFind3ds.dat'   using 1:8 with lines lt 3 lw 1 lc rgb "black"  title " {/Times-Italic l}_{max} = 10 ", \
     '3Dvlt_macro_noifthen_T_Cc0.40_lmax10_dl0.001_Op1_usingK0cFind1.dat'   using 1:8 with lines title " {/Times-Italic l}_{max} = 10 ", \
     '3Dvlt_macro_noifthen_T_Cc0.40_lmax10_dl0.001_Op1_usingK0cFind3dl.dat' using 1:8 with lines title " {/Times-Italic l}_{max} = 10 ", \
     '3Dvlt_macro_noifthen_T_Cc0.40_lmax10_dl0.001_Op1_usingK0cFind3ds.dat' using 1:8 with lines title " {/Times-Italic l}_{max} = 10 ", \
     '3Dvlt_macro_noifthen_T_Cc1.10_lmax10_dl0.001_Op1_usingK0cFind1.dat'   using 1:8 with lines title " {/Times-Italic l}_{max} = 10 ", \
     '3Dvlt_macro_noifthen_T_Cc1.10_lmax10_dl0.001_Op1_usingK0cFind3dl.dat' using 1:8 with lines title " {/Times-Italic l}_{max} = 10 ", \
     '3Dvlt_macro_noifthen_T_Cc1.10_lmax10_dl0.001_Op1_usingK0cFind3ds.dat' using 1:8 with lines title " {/Times-Italic l}_{max} = 10 ", \
     '3Dvlt_macro_noKth_T_Cc0.40_lmax10_dl0.001_Op1_usingK0cFind1.dat'     using 1:8 with lines lt 1 lw 3 lc rgb "yellow" title " {/Times-Italic l}_{max} = 10 ", \
     '3Dvlt_macro_noKth_T_Cc0.40_lmax10_dl0.001_Op1_usingK0cFind3dl.dat'   using 1:8 with lines lt 1 lw 3 lc rgb "yellow" title " {/Times-Italic l}_{max} = 10 ", \
     '3Dvlt_macro_noKth_T_Cc0.40_lmax10_dl0.001_Op1_usingK0cFind3ds.dat'   using 1:8 with lines lt 1 lw 3 lc rgb "yellow" title " {/Times-Italic l}_{max} = 10 ", \
     '3Dvlt_macro_noKth_T_Cc1.10_lmax10_dl0.001_Op1_usingK0cFind1.dat'     using 1:8 with lines lt 1 lw 3 lc rgb "yellow" title " {/Times-Italic l}_{max} = 10 ", \
     '3Dvlt_macro_noKth_T_Cc1.10_lmax10_dl0.001_Op1_usingK0cFind3dl.dat'   using 1:8 with lines lt 1 lw 3 lc rgb "yellow" title " {/Times-Italic l}_{max} = 10 ", \
     '3Dvlt_macro_noKth_T_Cc1.10_lmax10_dl0.001_Op1_usingK0cFind3ds.dat'   using 1:8 with lines lt 1 lw 3 lc rgb "yellow" title " {/Times-Italic l}_{max} = 10 "
     

