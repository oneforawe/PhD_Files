#filename: plot_3Dvlt_constT_Cc1.11_lmax10_dl0.001_Op1_acVsl.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_3Dvlt_constT_Cc1.11_lmax10_dl0.001_Op1_acVsl.ps'
#set size 1.0,1.1

set title "^4He Superfluid Vortex Loop Effective Core Diameter Fraction {/Times-Italic a}_c/{/Times-Italic a} vs Vortex Loop Effective Diameter {/Times-Italic a} \n [[this plot is not accurate... it does not correspond to how things are actually treated in the program]] \n [[the program should be re-written to explicitly calculate ac and a]]"

set ylabel "{/Times-Italic a}_c/{/Times-Italic a} (unitless)"
set xlabel "{/Times-Italic a} (unitless)"
#set xrange [0:1]
####set yrange [0:10]
#set logscale y
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
#set rmargin 35

plot '3Dvlt_constT_1_Cc1.11_lmax10_dl0.001_Op1.dat'   using (exp($3)):($4**0.6*exp(1.10501)) with linespoints title " lmax = 10 (correct as long as {/Times-Italic K}<1) ", \
     '3Dvlt_constT_0.9_Cc1.11_lmax10_dl0.001_Op1.dat' using (exp($3)):($4**0.6*exp(1.10501)) with linespoints title " lmax = 10 (correct as long as {/Times-Italic K}<1) ", \
     '3Dvlt_constT_0.8_Cc1.11_lmax10_dl0.001_Op1.dat' using (exp($3)):($4**0.6*exp(1.10501)) with linespoints title " lmax = 10 (correct as long as {/Times-Italic K}<1) ", \
     '3Dvlt_constT_0.7_Cc1.11_lmax10_dl0.001_Op1.dat' using (exp($3)):($4**0.6*exp(1.10501)) with linespoints title " lmax = 10 (correct as long as {/Times-Italic K}<1) ", \
     '3Dvlt_constT_0.6_Cc1.11_lmax10_dl0.001_Op1.dat' using (exp($3)):($4**0.6*exp(1.10501)) with linespoints title " lmax = 10 (correct as long as {/Times-Italic K}<1) ", \
     '3Dvlt_constT_0.5_Cc1.11_lmax10_dl0.001_Op1.dat' using (exp($3)):($4**0.6*exp(1.10501)) with linespoints title " lmax = 10 (correct as long as {/Times-Italic K}<1) ", \
     '3Dvlt_constT_0.4_Cc1.11_lmax10_dl0.001_Op1.dat' using (exp($3)):($4**0.6*exp(1.10501)) with linespoints title " lmax = 10 (correct as long as {/Times-Italic K}<1) ", \
     '3Dvlt_constT_0.3_Cc1.11_lmax10_dl0.001_Op1.dat' using (exp($3)):($4**0.6*exp(1.10501)) with linespoints title " lmax = 10 (correct as long as {/Times-Italic K}<1) ", \
     '3Dvlt_constT_0.2_Cc1.11_lmax10_dl0.001_Op1.dat' using (exp($3)):($4**0.6*exp(1.10501)) with linespoints title " lmax = 10 (correct as long as {/Times-Italic K}<1) ", \
     '3Dvlt_constT_0.1_Cc1.11_lmax10_dl0.001_Op1.dat' using (exp($3)):($4**0.6*exp(1.10501)) with linespoints title " lmax = 10 (correct as long as {/Times-Italic K}<1) "

#     '3Dvlt_constT_1_Cc1.11_lmax10_dl0.001_Op1.dat' using 3:4 with linespoints title " ({/Times-Italic K}) ", \
#     '3Dvlt_constT_0.9_Cc1.11_lmax10_dl0.001_Op1.dat' using 3:4 with linespoints title " ({/Times-Italic K}) ", \
#     '3Dvlt_constT_0.8_Cc1.11_lmax10_dl0.001_Op1.dat' using 3:4 with linespoints title " ({/Times-Italic K}) ", \
#     '3Dvlt_constT_0.7_Cc1.11_lmax10_dl0.001_Op1.dat' using 3:4 with linespoints title " ({/Times-Italic K}) ", \
#     '3Dvlt_constT_0.6_Cc1.11_lmax10_dl0.001_Op1.dat' using 3:4 with linespoints title " ({/Times-Italic K}) ", \
#     '3Dvlt_constT_0.5_Cc1.11_lmax10_dl0.001_Op1.dat' using 3:4 with linespoints title " ({/Times-Italic K}) ", \
#     '3Dvlt_constT_0.4_Cc1.11_lmax10_dl0.001_Op1.dat' using 3:4 with linespoints title " ({/Times-Italic K}) ", \
#     '3Dvlt_constT_0.3_Cc1.11_lmax10_dl0.001_Op1.dat' using 3:4 with linespoints title " ({/Times-Italic K}) ", \
#     '3Dvlt_constT_0.2_Cc1.11_lmax10_dl0.001_Op1.dat' using 3:4 with linespoints title " ({/Times-Italic K}) ", \
#     '3Dvlt_constT_0.1_Cc1.11_lmax10_dl0.001_Op1.dat' using 3:4 with linespoints title " ({/Times-Italic K}) "
