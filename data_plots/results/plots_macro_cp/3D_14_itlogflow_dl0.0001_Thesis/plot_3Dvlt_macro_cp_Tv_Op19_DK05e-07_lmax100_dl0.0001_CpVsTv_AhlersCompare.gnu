#filename: plot_3Dvlt_macro_cp_Tv_Op19_DK05e-07_lmax100_dl0.0001_CpVsTv_AhlersCompare.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_3Dvlt_macro_cp_Tv_Op19_DK05e-07_lmax100_dl0.0001_CpVsTv_AhlersCompare.ps'

set title "^4He Molar Specific Heat Capacity at Constant Pressure {/Times-Italic c_P} \n vs Temperature Variable {/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l}"

set ylabel "{/Times-Italic c_P} (J mol^{-1} K^{-1})"
set xlabel "{/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l}{/Times-Italic (P)} (unitless)"
set logscale x

## set xrange [0.1:1e-6]
## set yrange [0:240]
####### set xrange [0.1:5e-7]
set xrange [1:5e-7]
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

plot '3Dvlt_macro_cp_Tv_Op19_P00.050_DK05e-07_lmax100_dl0.0001.dat' using 4:(-1-$14) with linespoints notitle, \
     '3Dvlt_macro_cp_Tv_Op19_P00.050_DK05e-07_lmax100_dl0.0001.dat' using 4:14 with linespoints title " {/Times-Italic P} = {/Times-Italic P}_{sv} = 0.050 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.172 K) ", \
     '3Dvlt_macro_cp_Tv_Op19_P01.646_DK05e-07_lmax100_dl0.0001.dat' using 4:14 with linespoints title " {/Times-Italic P} = 1.644 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.157 K) ", \
     '3Dvlt_macro_cp_Tv_Op19_P07.328_DK05e-07_lmax100_dl0.0001.dat' using 4:14 with linespoints title " {/Times-Italic P} = 7.328 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.095 K) ", \
     '3Dvlt_macro_cp_Tv_Op19_P15.031_DK05e-07_lmax100_dl0.0001.dat' using 4:14 with linespoints title " {/Times-Italic P} = 15.031 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 1.998 K) ", \
     '3Dvlt_macro_cp_Tv_Op19_P18.180_DK05e-07_lmax100_dl0.0001.dat' using 4:14 with linespoints title " {/Times-Italic P} = 18.180 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 1.954 K) ", \
     '3Dvlt_macro_cp_Tv_Op19_P25.868_DK05e-07_lmax100_dl0.0001.dat' using 4:14 with linespoints title " {/Times-Italic P} = 25.868 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 1.836 K) ", \
     -1397.31440791919*x**0.0150650568395361 + 1371.98625184181 ls 0 title " Fits ", \
     -1406.29005656644*x**0.0150650568395361 + 1381.03734232903 ls 0 notitle, \
     -1400.57589816968*x**0.0150650568395361 + 1377.12888985157 ls 0 notitle, \
     -1327.60605121375*x**0.0150650568395361 + 1308.41112140293 ls 0 notitle, \
     -1285.62078993663*x**0.0150650568395361 + 1268.63354185102 ls 0 notitle, \
     -1173.01667705812*x**0.0150650568395361 + 1163.60322249204 ls 0 notitle

#     -18.2455538115244*log(x) + -12.6811348796468 lt 0 lc rgb "red" title " Fits ", \
#     -18.3508003559673*log(x) + -12.4258422434769 lt 0 lc rgb "red" notitle, \
#     -18.272480072844*log(x) + -10.6460182453261  lt 0 lc rgb "red" notitle, \
#     -17.3038654985909*log(x) + -6.90815905803329 lt 0 lc rgb "red" notitle, \
#     -16.7834957276909*log(x) + -5.34152200759794 lt 0 lc rgb "red" notitle, \
#     -15.9547466616406*log(x) + -2.00726292281082 lt 0 lc rgb "red" notitle, \
#     -15.3128990132488*log(x) + 1.23310037226038  lt 0 lc rgb "red" notitle



#     '3Dvlt_macro_cp_Tv_Op19_P22.533_DK05e-07_lmax100_dl0.0001.dat' using 4:14 with linespoints title " {/Times-Italic P} = 22.533 bar ({/Times-Italic T}_{/Symbol-Oblique l} = ..... K) ", \
#     -1223.36373708572*x**0.0150650568395361 + 1210.10815577063 ls 0 notitle, \
