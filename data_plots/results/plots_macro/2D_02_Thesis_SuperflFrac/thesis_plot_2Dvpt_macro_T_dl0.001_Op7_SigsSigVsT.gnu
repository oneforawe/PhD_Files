#filename: thesis_plot_2Dvpt_macro_T_dl0.001_Op7_SigsSigVsT.gnu
reset
set terminal postscript color eps enhanced
set output 'thesis_plot_2Dvpt_macro_T_dl0.001_Op7_SigsSigVsT.ps'
set size 0.85,0.85

set title "^4He Film Superfluid Fraction {/Symbol-Oblique s}_s / {/Symbol-Oblique s} vs Temperature Fraction {/Times-Italic T} / {/Times-Italic T}_{KT}"

set ylabel "{/Symbol-Oblique s}_s / {/Symbol-Oblique s} (unitless)"
set xlabel "{/Times-Italic T} / {/Times-Italic T}_{KT} (unitless)"
set xrange [0:2]
set yrange [0:1.1]
#set logscale y
#set grid
#set key inside left bottom box width -2 height 0.5 spacing 1.2
#set key at graph 0.325,0.42 box width -2 height 0.5 spacing 1.2 title " ^{}Vortex Pair Theory "
set key at graph 0.325,0.53 box width -2 height 0.5 spacing 1.2 title " ^{}Vortex Pair Theory "
#set rmargin 35

#set pointsize 0.5

# lmax=20 (~1 cm), K/K0 @ Tkt
set ytics add ("" 0.869702)

# lmax=83 (~10^26 m, size of universe), K/K0 @ Tkt
set ytics add ("" 0.857343)

# lmax=100, K/K0 @ Tkt
#set ytics add ("" 0.856744)

# lmax= infty, K/K0 @ Tkt
set ytics add ("" (2.0/pi)/0.747853)



plot '2Dvpt_macro_T_lmax0_dl-nan_Op7.dat'    using 1:5 with lines title " {/Times-Italic l}_{max} = 0 ", \
     '2Dvpt_macro_T_lmax1_dl0.001_Op7.dat'   using 1:5 with lines title " {/Times-Italic l}_{max} = 1 ", \
     '2Dvpt_macro_T_lmax2_dl0.001_Op7.dat'   using 1:5 with lines title " {/Times-Italic l}_{max} = 2 ", \
     '2Dvpt_macro_T_lmax4_dl0.001_Op7.dat'   using 1:5 with lines title " {/Times-Italic l}_{max} = 4 ", \
     '2Dvpt_macro_T_lmax10_dl0.001_Op7.dat'  using 1:5 with lines title " {/Times-Italic l}_{max} = 10 ", \
     '2Dvpt_macro_T_lmax16_dl0.001_Op7.dat'  using 1:5 with lines title " {/Times-Italic l}_{max} = 16 ", \
     '2Dvpt_macro_T_lmax23_dl0.001_Op7.dat'  using 1:5 with lines title " {/Times-Italic l}_{max} = 23 ", \
     '2Dvpt_macro_T_lmax83_dl0.001_Op7.dat'  using 1:5 with lines title " {/Times-Italic l}_{max} = 83 "
#     0.856744 ls 0 notitle

#lc rgb "black" 

#     '2Dvpt_macro_T_lmax25_dl0.001_Op7.dat'  using 1:5 with lines title " {/Times-Italic l}_{max} = 25 ", \
#     '2Dvpt_macro_T_lmax22_dl0.001_Op7.dat'  using 1:5 with lines title " {/Times-Italic l}_{max} = 22 ", \
#     '2Dvpt_macro_T_lmax10_dl0.001_Op7.dat'  using 1:5 with lines title " {/Times-Italic l}_{max} = 20 ", \
#     '2Dvpt_macro_T_lmax30_dl0.001_Op7.dat'  using 1:5 with lines title " {/Times-Italic l}_{max} = 30 ", \
#     '2Dvpt_macro_T_lmax37_dl0.001_Op7.dat'  using 1:5 with lines title " {/Times-Italic l}_{max} = 37 ", \
#     '2Dvpt_macro_T_lmax50_dl0.001_Op7.dat'  using 1:5 with lines title " {/Times-Italic l}_{max} = 50 ", \
#     '2Dvpt_macro_T_lmax90_dl0.001_Op7.dat'  using 1:5 with lines title " {/Times-Italic l}_{max} = 90 ", \
#     '2Dvpt_macro_T_lmax100_dl0.001_Op7.dat' using 1:5 with lines title " {/Times-Italic l}_{max} = 100 ", \
#     '2Dvpt_macro_T_lmax175_dl0.001_Op7.dat' using 1:5 with lines title " {/Times-Italic l}_{max} = 175 ", \


#     '2Dvpt_macro_T_lmax3_dl0.001_Op7.dat'   using 1:5 with lines title " {/Times-Italic l}_{max} = 3 ", \
#     '2Dvpt_macro_T_lmax4_dl0.001_Op7.dat'   using 1:5 with lines title " {/Times-Italic l}_{max} = 4 ", \
#     '2Dvpt_macro_T_lmax6_dl0.001_Op7.dat'   using 1:5 with lines title " {/Times-Italic l}_{max} = 6 ", \
#     '2Dvpt_macro_T_lmax7_dl0.001_Op7.dat'   using 1:5 with lines title " {/Times-Italic l}_{max} = 7 ", \
#     '2Dvpt_macro_T_lmax8_dl0.001_Op7.dat'   using 1:5 with lines title " {/Times-Italic l}_{max} = 8 ", \
#     '2Dvpt_macro_T_lmax9_dl0.001_Op7.dat'   using 1:5 with lines title " {/Times-Italic l}_{max} = 9 ", \
#     '2Dvpt_macro_T_lmax9_dl0.001_Op7_K0c0.74785238.dat'   using 1:5 with lines title " l_{max} = 9* ", \
