#filename: plot_3Dvlt_macro_cp_Tv_Op42_Kth4.147_DK05e-07_lmax100_dl1e-05_CpVsTv_AhlersCompare.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_3Dvlt_macro_cp_Tv_Op42_Kth4.147_DK05e-07_lmax100_dl1e-05_CpVsTv_AhlersCompare.ps'

set title "^4He Molar Specific Heat Capacity at Constant Pressure {/Times-Italic c_P} \n vs Temperature Variable {/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l}"

set ylabel "{/Times-Italic c_P} (J mol^{-1} K^{-1})"
set xlabel "{/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l}{/Times-Italic (P)} (unitless)"
set logscale x

## set xrange [0.1:1e-6]
## set yrange [0:240]
set xrange [0.1:5e-7]
set yrange [0:300]
#set xrange [1:1e-8]
#set yrange [0:300]
##set xrange [1e-2:1e-8]
set mxtics 10

set logscale x
### set xrange [1:1e-150]
### set yrange [0:1450]
#set yrange [0:450]

set grid
#set key inside left top box width -8 spacing 1.2 title "Calculated from \n {/CM-Typewriter vlt\\_HeatCap\\_P.c}"
set key inside left top box width -8 height 0.5 spacing 1.3 title " VLT calculations "
set rmargin 4

plot '3Dvlt_macro_cp_Tv_Op42_Kth4.147_P00.050_DK05e-07_lmax100_dl1e-05.dat' using 4:(-1-$14) with linespoints notitle, \
     '3Dvlt_macro_cp_Tv_Op42_Kth4.147_P00.050_DK05e-07_lmax100_dl1e-05.dat' using 4:14 with linespoints title " {/Times-Italic P} = {/Times-Italic P}_{sv} = 0.050 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.172 K) ", \
     '3Dvlt_macro_cp_Tv_Op42_Kth4.147_P01.646_DK05e-07_lmax100_dl1e-05.dat' using 4:14 with linespoints title " {/Times-Italic P} = 1.644 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.157 K) ", \
     '3Dvlt_macro_cp_Tv_Op42_Kth4.147_P07.328_DK05e-07_lmax100_dl1e-05.dat' using 4:14 with linespoints title " {/Times-Italic P} = 7.328 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.095 K) ", \
     '3Dvlt_macro_cp_Tv_Op42_Kth4.147_P15.031_DK05e-07_lmax100_dl1e-05.dat' using 4:14 with linespoints title " {/Times-Italic P} = 15.031 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 1.998 K) ", \
     '3Dvlt_macro_cp_Tv_Op42_Kth4.147_P18.180_DK05e-07_lmax100_dl1e-05.dat' using 4:14 with linespoints title " {/Times-Italic P} = 18.180 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 1.954 K) ", \
     '3Dvlt_macro_cp_Tv_Op42_Kth4.147_P25.868_DK05e-07_lmax100_dl1e-05.dat' using 4:14 with linespoints title " {/Times-Italic P} = 25.868 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 1.836 K) ", \
     -1426.65983839276*x**0.0150650568395361 + 1401.04566791893 ls 0 title " Fits ", \
     -1435.92303060133*x**0.0150650568395361 + 1410.36520753 ls 0 notitle, \
     -1429.8538333734*x**0.0150650568395361 + 1406.08470559295 ls 0 notitle, \
     -1354.47581410578*x**0.0150650568395361 + 1334.99373819771 ls 0 notitle, \
     -1312.66694007668*x**0.0150650568395361 + 1295.24544111872 ls 0 notitle, \
     -1198.99819642038*x**0.0150650568395361 + 1188.95719626103 ls 0 notitle


#     '3Dvlt_macro_cp_Tv_Op42_Kth4.147_P22.533_DK05e-07_lmax100_dl1e-05.dat' using 4:14 with linespoints title " {/Times-Italic P} = 22.533 bar ({/Times-Italic T}_{/Symbol-Oblique l} = ..... K) ", \
#     -1248.23211337185*x**0.0150650568395361 + 1234.63543285753 ls 0 notitle, \
