#filename: plot_3Dvlt_macro_cp_DK0_5e-07_i47_PartlyCorrected_CpVsTv_AhlersCompare.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_3Dvlt_macro_cp_DK0_5e-07_i47_PartlyCorrected_CpVsTv_AhlersCompare.ps'

set title "^4He Molar Specific Heat Capacity at Constant Pressure {/Times-Italic c_P} \n vs Temperature Variable {/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l}"

set ylabel "{/Times-Italic c_P} (J mol^{-1} K^{-1})"
set xlabel "{/Symbol-Oblique t} = {/Times-Italic 1 - T / T}_{/Symbol-Oblique l}{/Times-Italic (P)} (unitless)"
set logscale x
set xrange [1:1e-8]
#set xrange [1e-2:1e-8]
set mxtics 10
set yrange [0:300]
set grid
#set key inside left top box width -8 spacing 1.2 title "Calculated from \n {/CM-Typewriter vlt\\_HeatCap\\_P.c}"
set key inside left top box width -8 height 0.5 spacing 1.3 title " VLT calculations "
set rmargin 4

plot '3Dvlt_macro_cp_P00.050_DK05e-07_i47_PartlyCorrected.dat' using 2:(-1-$9) with linespoints notitle, \
     '3Dvlt_macro_cp_P00.050_DK05e-07_i47_PartlyCorrected.dat' using 2:9 with linespoints title " {/Times-Italic P} = {/Times-Italic P}_{sv} = 0.050 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.172 K) ", \
     '3Dvlt_macro_cp_P01.646_DK05e-07_i47_PartlyCorrected.dat' using 2:9 with linespoints title " {/Times-Italic P} = 1.644 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.157 K) ", \
     '3Dvlt_macro_cp_P07.328_DK05e-07_i47_PartlyCorrected.dat' using 2:9 with linespoints title " {/Times-Italic P} = 7.328 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.095 K) ", \
     '3Dvlt_macro_cp_P15.031_DK05e-07_i47_PartlyCorrected.dat' using 2:9 with linespoints title " {/Times-Italic P} = 15.031 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 1.998 K) ", \
     '3Dvlt_macro_cp_P18.180_DK05e-07_i47_PartlyCorrected.dat' using 2:9 with linespoints title " {/Times-Italic P} = 18.180 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 1.954 K) ", \
     '3Dvlt_macro_cp_P25.868_DK05e-07_i47_PartlyCorrected.dat' using 2:9 with linespoints title " {/Times-Italic P} = 25.868 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 1.836 K) ", \
     -1423.54077092243*x**0.0150649 + 1397.98033663787 ls 0 title " Fits ", \
     -1433.00550290483*x**0.0150649 + 1406.68877142855 ls 0 notitle, \
     -1426.45092191267*x**0.0150649 + 1399.40613082434 ls 0 notitle, \
     -1352.11655226025*x**0.0150649 + 1325.4633849003  ls 0 notitle, \
     -1310.39280789169*x**0.0150649 + 1283.75091991116 ls 0 notitle, \
     -1194.47697925562*x**0.0150649 + 1169.28538526082 ls 0 notitle

#     '3Dvlt_macro_cp_P22.533_DK05e-07_i47_PartlyCorrected.dat' using 2:9 with linespoints title " {/Times-Italic P} = 22.533 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 1.889 K) ", \
#     -1246.32596847999*x**0.0150649 + 1220.05569788508 notitle, \
