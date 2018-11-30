#filename: thesis_plot_2Dvpt_macro_T_lmax_dl0.001_Op7_SigsSigVsT.gnu
reset
set terminal postscript color eps enhanced

set output 'thesis_plot_2Dvpt_macro_T_lmax5_dl0.001_Op7_SigsSigVsT_v2.ps'
set size 0.85,0.88
#set xrange [0:1.07]
set xrange [0:1.2]

#set output 'thesis_plot_2Dvpt_macro_T_lmax4and5_dl0.001_Op7_SigsSigVsT.ps'
#set size 0.85,0.88
#set xrange [0:1.2]

#set output 'thesis_plot_2Dvpt_macro_T_lmax5_dl0.001_Op7_SigsSigVsT.ps'
#set size 0.85,0.88
#set xrange [0:1.2]

#set output 'thesis_plot_2Dvpt_macro_T_lmax4_dl0.001_Op7_SigsSigVsT.ps'
#set size 0.84,0.88
#set xrange [0:0.8]

#set output 'thesis_plot_2Dvpt_macro_T_lmax3.5_dl0.001_Op7_SigsSigVsT.ps'
#set size 0.84,0.88
#set xrange [0:1.2]

set title "^4He Film Superfluid Fraction {/Symbol-Oblique s}_s / {/Symbol-Oblique s} vs Temperature Fraction {/Times-Italic T} / {/Times-Italic T}_{KT}"

set ylabel "{/Symbol-Oblique s}_s / {/Symbol-Oblique s} (unitless)"
set xlabel "{/Times-Italic T} / {/Times-Italic T}_{KT} (unitless)"
set yrange [0:1.1]
#set logscale y
#set grid
#set key inside left bottom box width -2 height 0.5 spacing 1.2
#set key at graph 0.325,0.17 box width -2 height 0.5 spacing 1.2 title " ^{}Vortex Pair Theory "
set key at graph 0.325,0.2 box width -6 height 0.5 spacing 1.2 title " ^{}Vortex Pair Theory, \n ^{}modified Villain model "
#set rmargin 35

#set pointsize 0.5

plot '2Dvpt_macro_T_lmax5_dl0.001_Op7.dat'             using 1:5 with lines lc rgb "blue" lt 2 lw 2 title " {/Times-Italic l}_{max} = 5 ", \
     '2Dvpt_macro_T_lmax5_dl0.001_Op7_Villain0.8.dat'  using 1:5 with lines lc rgb "green" lt 2 lw 2 title " {/Times-Italic l}_{max} = 5 mVm0.8 ", \
     '2Dvpt_macro_T_lmax5_dl0.001_Op7_Villain0.75.dat' using 1:5 with lines lc rgb "green" lt 2 lw 2 title " {/Times-Italic l}_{max} = 5 mVm0.75 ", \
     '2Dvpt_macro_T_lmax5_dl0.001_Op7_VillainMod.dat'  using 1:5 with lines lc rgb "green" lt 2 lw 2 title " {/Times-Italic l}_{max} = 5 mVm0.5 ", \
     '2Dvpt_macro_T_lmax5_dl0.001_Op7_Villain0.25.dat' using 1:5 with lines lc rgb "green" lt 2 lw 2 title " {/Times-Italic l}_{max} = 5 mVm0.25 ", \
     '2Dvpt_macro_T_lmax5_dl0.001_Op7_Villain0.1.dat'  using 1:5 with lines lc rgb "green" lt 2 lw 2 title " {/Times-Italic l}_{max} = 5 mVm0.1 "
#     '2Dvpt_macro_T_lmax4_dl0.001_Op7_Villain0.8.dat' using 1:5 with lines lc rgb "green" lt 2 lw 2 title " {/Times-Italic l}_{max} = 4 mVm0.8 "
#     '2Dvpt_macro_T_lmax3.5_dl0.001_Op7_Villain0.75.dat' using 1:5 with lines lc rgb "green" lt 2 lw 2 title " {/Times-Italic l}_{max} = 3.5 mVm0.75 ", \
#     '2Dvpt_macro_T_lmax4_dl0.001_Op7_VillainMod.dat' using 1:5 with lines lc rgb "green" lt 2 lw 2 title " {/Times-Italic l}_{max} = 4 mVm0.5", \

#     The plot above "lmax3.5" uses the modified Villain model ("mVm" 0.75* )
#     The plot above "lmax4"   uses the modified Villain model ("mVm" 0.5* )


#gs -sDEVICE=png16 -sOutputFile=outfile.png -r300 infile.ps
#gs -sDEVICE=png16 -sOutputFile=thesis_plot_2Dvpt_macro_T_lmax5_dl0.001_Op7_SigsSigVsT.png -r300 thesis_plot_2Dvpt_macro_T_lmax5_dl0.001_Op7_SigsSigVsT.ps
#gs -sDEVICE=png16 -sOutputFile=thesis_plot_2Dvpt_macro_T_lmax5_dl0.001_Op7_SigsSigVsT_b.png -r900 thesis_plot_2Dvpt_macro_T_lmax5_dl0.001_Op7_SigsSigVsT.ps

#gs -sDEVICE=png16 -sOutputFile=thesis_plot_2Dvpt_macro_T_lmax4_dl0.001_Op7_SigsSigVsT.png -r900 thesis_plot_2Dvpt_macro_T_lmax4_dl0.001_Op7_SigsSigVsT.ps
#gs -sDEVICE=png16 -sOutputFile=thesis_plot_2Dvpt_macro_T_lmax3.5_dl0.001_Op7_SigsSigVsT.png -r900 thesis_plot_2Dvpt_macro_T_lmax3.5_dl0.001_Op7_SigsSigVsT.ps
