#filename: plot_3Dvlt_macro_cp_Tv_Op19_DK05e-07_lmax100_dl0.0001_CpVsTv_TermFrac.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_3Dvlt_macro_cp_Tv_Op19_DK05e-07_lmax100_dl0.0001_CpVsTv_TermFrac.ps'

set title "^4He Molar Specific Heat Capacity Terms {/Times-Italic c_P} = cp01+cp02+... \n vs Temperature Variable {/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l}"

set ylabel "{/Times-Italic c_P} terms (J mol^{-1} K^{-1})"
set xlabel "{/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l}{/Times-Italic (P)} (unitless)"
set logscale x
set xrange [1:1e-8]
#set xrange [1e-2:1e-8]
set mxtics 10
#set yrange [0:300]
#set yrange [-30:300]
#set yrange [-0.1:1.2]
## set yrange [-1:3]
#set yrange [-1:20]
set grid
#set key inside left top box width -8 spacing 1.2 title "Calculated from \n {/CM-Typewriter vlt\\_HeatCap\\_P.c}"
set key outside right top box width -10 height 0.5 spacing 1.3 title " VLT calculations "
set rmargin 20

plot '3Dvlt_macro_cp_Tv_Op19_P00.050_DK05e-07_lmax100_dl0.0001.dat' using 4:($19/$14)  with linespoints title " 00.050 bar; cp05/cp ", \
     '3Dvlt_macro_cp_Tv_Op19_P01.646_DK05e-07_lmax100_dl0.0001.dat' using 4:($19/$14)  with linespoints title " 01.646 bar; cp05/cp ", \
     '3Dvlt_macro_cp_Tv_Op19_P07.328_DK05e-07_lmax100_dl0.0001.dat' using 4:($19/$14)  with linespoints title " 07.328 bar; cp05/cp ", \
     '3Dvlt_macro_cp_Tv_Op19_P15.031_DK05e-07_lmax100_dl0.0001.dat' using 4:($19/$14)  with linespoints title " 15.031 bar; cp05/cp ", \
     '3Dvlt_macro_cp_Tv_Op19_P18.180_DK05e-07_lmax100_dl0.0001.dat' using 4:($19/$14)  with linespoints title " 18.180 bar; cp05/cp ", \
     '3Dvlt_macro_cp_Tv_Op19_P22.533_DK05e-07_lmax100_dl0.0001.dat' using 4:($19/$14)  with linespoints title " 22.533 bar; cp05/cp ", \
     '3Dvlt_macro_cp_Tv_Op19_P25.868_DK05e-07_lmax100_dl0.0001.dat' using 4:($19/$14)  with linespoints title " 25.868 bar; cp05/cp "
