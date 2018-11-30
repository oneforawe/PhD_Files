#filename: plot_3Dvlt_macro_cp_DK0_5e-07_i47_CpVsTv_AhlersCompare.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_3Dvlt_macro_cp_DK0_5e-07_i47_CpVsTv_AhlersCompare.ps'

set title "^4He Molar Specific Heat Capacity at Constant Pressure {/Times-Italic c_P} \n vs Temperature Variable {/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l}"

set ylabel "{/Times-Italic c_P} (J mol^{-1} K^{-1})"
set xlabel "{/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l}{/Times-Italic (P)} (unitless)"
set logscale x
set xrange [0.1:1e-6]
#set xrange [1:1e-8]
##set xrange [1e-2:1e-8]
set mxtics 10
#set yrange [0:300]
set yrange [0:240]
set grid
#set key inside left top box width -8 spacing 1.2 title "Calculated from \n {/CM-Typewriter vlt\\_HeatCap\\_P.c}"
set key inside left top box width -8 height 0.5 spacing 1.3 title " VLT calculations "
set rmargin 4

plot '3Dvlt_macro_cp_P00.050_DK05e-07_i47.dat' using 2:(-1-$9) with linespoints notitle, \
     '3Dvlt_macro_cp_P00.050_DK05e-07_i47.dat' using 2:9 with linespoints title " {/Times-Italic P} = {/Times-Italic P}_{sv} = 0.050 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.172 K) ", \
     '3Dvlt_macro_cp_P01.646_DK05e-07_i47.dat' using 2:9 with linespoints title " {/Times-Italic P} = 1.644 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.157 K) ", \
     '3Dvlt_macro_cp_P07.328_DK05e-07_i47.dat' using 2:9 with linespoints title " {/Times-Italic P} = 7.328 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.095 K) ", \
     '3Dvlt_macro_cp_P15.031_DK05e-07_i47.dat' using 2:9 with linespoints title " {/Times-Italic P} = 15.031 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 1.998 K) ", \
     '3Dvlt_macro_cp_P18.180_DK05e-07_i47.dat' using 2:9 with linespoints title " {/Times-Italic P} = 18.180 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 1.954 K) ", \
     '3Dvlt_macro_cp_P25.868_DK05e-07_i47.dat' using 2:9 with linespoints title " {/Times-Italic P} = 25.868 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 1.836 K) ", \
     -1422.52535699295*x**0.0150650568 + 1391.91115973955 ls 0 title " Fits ", \
     -1432.06707211571*x**0.0150650568 + 1400.58269721735 ls 0 notitle, \
     -1424.77737006795*x**0.0150650568 + 1388.75150340618 ls 0 notitle, \
     -1349.27760664466*x**0.0150650568 + 1308.7485893577  ls 0 notitle, \
     -1307.44033540544*x**0.0150650568 + 1266.70091179144 ls 0 notitle, \
     -1190.85136252856*x**0.0150650568 + 1149.5938497745  ls 0 notitle

#     -1243.28354610711*x**0.0150650568 + 1202.95368702857 ls 0 notitle, \

#     Incorrect (but oh-so-nice) fits using "20-term" cp
#     -1323.90265363016*x**0.0150650568 + 1300.70053643625 ls 0 title " Fits ", \
#     -1308.17644364677*x**0.0150650568 + 1284.81609410076 ls 0 notitle, \
#     -1224.31932856669*x**0.0150650568 + 1202.42212629271 ls 0 notitle, \
#     -1114.57503185628*x**0.0150650568 + 1094.63208887831 ls 0 notitle, \
#     -1072.82579355647*x**0.0150650568 + 1053.07811648916 ls 0 notitle, \
#     -953.782056115115*x**0.0150650568 + 936.060503397096 ls 0 notitle

#     Fits from prior to the addition of ten more terms in cp
#     -1423.54077092243*x**0.0150649 + 1397.98033663787 ls 0 title " Fits ", \
#     -1433.00550290483*x**0.0150649 + 1406.68877142855 ls 0 notitle, \
#     -1426.45092191267*x**0.0150649 + 1399.40613082434 ls 0 notitle, \
#     -1352.11655226025*x**0.0150649 + 1325.4633849003  ls 0 notitle, \
#     -1310.39280789169*x**0.0150649 + 1283.75091991116 ls 0 notitle, \
#     -1194.47697925562*x**0.0150649 + 1169.28538526082 ls 0 notitle
