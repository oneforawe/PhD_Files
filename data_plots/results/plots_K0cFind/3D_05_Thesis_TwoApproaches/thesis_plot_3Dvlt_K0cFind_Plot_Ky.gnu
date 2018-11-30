#filename: thesis_plot_3Dvlt_K0cFind_Plot_Ky.gnu
reset
set terminal postscript color eps enhanced
set output 'thesis_plot_3Dvlt_K0cFind_Plot_Ky.ps'
#set size 0.85,0.85
##set size 1.00,0.85
set size 0.85,1.20

#set samples 1000000

set title "^4He Fugacity {/Times-Italic y} vs Coupling Parameter {/Times-Italic K}"

set ylabel "{/Times-Italic y} (unitless)"
set xlabel "{/Times-Italic K} (unitless)"

#set xrange [6:8]
#set yrange [0:0.01]

#set xrange [30.5:31.5]
#set yrange [-0.01:0.01]

#set xrange [-0.1:35.5]
#set yrange [-0.1:0.3]

set xrange [-0.1:2]
set yrange [-0.03:0.27]
#set xrange [-0.1:9]
#set yrange [-0.03:0.1]
#set xrange [-0.1:8]
#set yrange [-0.01:0.1]

#set xrange [-0.1:5.5]
#set yrange [-0.01:0.3]
#set xrange [-0.1:1]
#set yrange [-0.01:0.1]

##set xrange [0.25:0.5]
##set yrange [0.04:0.08]
###set xrange [0.2:1.5]
###set yrange [0.01:0.15]
#set xrange [0.387508189712343-3e-4:0.387508189712343+3e-4]
#set yrange [0.06242100545760196-5e-5:0.06242100545760196+5e-5]

##set xrange [0.387508189712343-3e-5:0.387508189712343+3e-5]
##set yrange [0.06242100545760196-5e-6:0.06242100545760196+5e-6]

###set xrange [0.387508189712343-3e-6:0.387508189712343+3e-6]
###set yrange [0.06242100545760196-5e-7:0.06242100545760196+5e-7]
####set xrange [0.387508189712343-1.1e-6:0.387508189712343+1.1e-6]
####set yrange [0.06242100545760196-4e-7:0.06242100545760196+4e-7]

#set xrange [0.387508189712343-3e-7:0.387508189712343+3e-7]
#set yrange [0.06242100545760196-5e-8:0.06242100545760196+5e-8]

#set xrange [0.387508189712343-3e-14:0.387508189712343+3e-14]
#set yrange [0.06242100545760196-5e-15:0.06242100545760196+5e-15]
#set xrange [0.387508189712343-3e-15:0.387508189712343+3e-15]
#set yrange [0.06242100545760196-5e-16:0.06242100545760196+5e-16]
#set xrange [0.387508189712343-3e-16:0.387508189712343+3e-16]
#set yrange [0.06242100545760196-5e-17:0.06242100545760196+5e-17]

##set xrange [3.4:5.2]
##set yrange [0.003:0.012]

set grid
#####set key at graph 0.85,-0.21 box width -5 height 0.5 spacing 1.6 title " Finding {/Times-Italic K}@^{{} {/Times-Italic l}0^{^{}}}_{c  }, given a {/Times-Italic C}_c^ "
##set key at graph 0.875,-0.21 box width -5 height 0.5 spacing 1.6 title " Finding {/Times-Italic K}@^{{} {/Times-Italic l}0^{^{}}}_{c  }, given a {/Times-Italic C}_c^ "
#set key outside below center box width -5 height 0.5 spacing 1.2 title " Finding {/Times-Italic K}_{0c}, given a {/Times-Italic C}@_c^ "
#####set bmargin 9.7

##set rmargin 20
##set key outside right
set bmargin 20
set key outside below

set pointsize 0.5

PISQ = 9.86960440108935861883
PICU = 31.00627668029981620634
A3 = 4.0*PICU/3.0

plot '3Dvlt_K0cFind2_Plot_Cc0.10_try10.dat' using 2:3 with points title " '2try10' 0.10 ", \
     '3Dvlt_K0cFind2_Plot_Cc1.20_try10.dat' using 2:3 with points title " '2try10' 1.20 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.10_try9.dat' using 2:3 with points title " '2try9' 0.10 ", \
     '3Dvlt_K0cFind2_Plot_Cc1.20_try9.dat' using 2:3 with points title " '2try9' 1.20 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.10_try7.dat' using 2:3 with points title " '2try7' 0.10 ", \
     '3Dvlt_K0cFind2_Plot_Cc1.20_try7.dat' using 2:3 with points title " '2try7' 1.20 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.10_try13.dat' using 2:3 with points title " '2try13' 0.10 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.10_try12.dat' using 2:3 with points title " '2try12' 0.10 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.10_try11.dat' using 2:3 with points title " '2try11' 0.10 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.10_try6.dat' using 2:3 with points title " '2try6' 0.10 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.10_try5.dat' using 2:3 with points title " '2try5' 0.10 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.20_try5.dat' using 2:3 with points title " '2try5' 0.20 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.30_try5.dat' using 2:3 with points title " '2try5' 0.30 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.40_try5.dat' using 2:3 with points title " '2try5' 0.40 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.50_try5.dat' using 2:3 with points title " '2try5' 0.50 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.60_try5.dat' using 2:3 with points title " '2try5' 0.60 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.70_try5.dat' using 2:3 with points title " '2try5' 0.70 ", \
     '3Dvlt_K0cFind2_Plot_Cc1.20_try6.dat' using 2:3 with points title " '2try6' 1.20 ", \
     '3Dvlt_K0cFind2_Plot_Cc1.20_try5.dat' using 2:3 with points title " '2try5' 1.20 ", \
     '3Dvlt_K0cFind2_Plot_Cc1.10_try5.dat' using 2:3 with points title " '2try5' 1.10 ", \
     '3Dvlt_K0cFind2_Plot_Cc1.00_try5.dat' using 2:3 with points title " '2try5' 1.00 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.90_try5.dat' using 2:3 with points title " '2try5' 0.90 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.80_try5.dat' using 2:3 with points title " '2try5' 0.80 ", \
     '3Dvlt_K0cFind2_Plot_Cc1.20_try7.dat' using 2:3 with points title " '2try7' 0.10 ", \
     '3Dvlt_K0cFind2_Plot_Cc1.20_try2.dat' using 2:3 with points title " '2try2' 0.10 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.10_try4.dat' using 2:3 with points title " '2try4' 0.10 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.20_try4.dat' using 2:3 with points title " '2try4' 0.20 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.30_try4.dat' using 2:3 with points title " '2try4' 0.30 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.40_try4.dat' using 2:3 with points title " '2try4' 0.40 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.50_try4.dat' using 2:3 with points title " '2try4' 0.50 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.60_try4.dat' using 2:3 with points title " '2try4' 0.60 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.70_try4.dat' using 2:3 with points title " '2try4' 0.70 ", \
     '3Dvlt_K0cFind2_Plot_Cc1.20_try7.dat' using 2:3 with points title " '2try7' 1.20 ", \
     '3Dvlt_K0cFind2_Plot_Cc1.20_try4.dat' using 2:3 with points title " '2try4' 1.20 ", \
     '3Dvlt_K0cFind2_Plot_Cc1.10_try4.dat' using 2:3 with points title " '2try4' 1.10 ", \
     '3Dvlt_K0cFind2_Plot_Cc1.00_try4.dat' using 2:3 with points title " '2try4' 1.00 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.90_try4.dat' using 2:3 with points title " '2try4' 0.90 ", \
     '3Dvlt_K0cFind2_Plot_Cc0.80_try4.dat' using 2:3 with points title " '2try4' 0.80 ", \
     '3Dvlt_K0cFind_Plot_Cc0.10_last.dat' using 2:3 with points title " 0.10 ", \
     exp(-PISQ*0.10*x) lt 2 lc rgb "yellow" notitle, \
     exp(-PISQ*0.20*x) lt 2 lc rgb "yellow" notitle, \
     exp(-PISQ*0.30*x) lt 2 lc rgb "yellow" notitle, \
     exp(-PISQ*0.40*x) lt 2 lc rgb "blue" notitle, \
     exp(-PISQ*0.50*x) lt 2 notitle, \
     exp(-PISQ*0.60*x) lt 2 notitle, \
     exp(-PISQ*0.70*x) lt 2 notitle, \
     exp(-PISQ*0.80*x) lt 2 notitle, \
     exp(-PISQ*0.90*x) lt 2 notitle, \
     exp(-PISQ*1.00*x) lt 2 notitle, \
     exp(-PISQ*1.10*x) lt 2 lc rgb "red"  notitle, \
     exp(-PISQ*1.20*x) lt 2 lc rgb "yellow" notitle, \
     1.0/(A3*x) notitle ls 0, \
     "<echo 0.387508189712343 0.06242100545760196" with points lc rgb "black" title " Fixed pts ", \
     "<echo 0.387508189712343409 0.0624210054576019" with points lc rgb "black" title " Fixed pts ", \
     "<echo 4.146766416946759293 0.0058331356032128735" with points lc rgb "black" notitle




#     '3Dvlt_K0cFind2_Plot_Cc0.10_try3.dat' using 2:3 with points title " '2try3' 0.10 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.20_try3.dat' using 2:3 with points title " '2try3' 0.20 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.30_try3.dat' using 2:3 with points title " '2try3' 0.30 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.40_try3.dat' using 2:3 with points title " '2try3' 0.40 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.50_try3.dat' using 2:3 with points title " '2try3' 0.50 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.60_try3.dat' using 2:3 with points title " '2try3' 0.60 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.70_try3.dat' using 2:3 with points title " '2try3' 0.70 ", \
#     '3Dvlt_K0cFind2_Plot_Cc1.20_try3.dat' using 2:3 with points title " '2try3' 1.20 ", \
#     '3Dvlt_K0cFind2_Plot_Cc1.10_try3.dat' using 2:3 with points title " '2try3' 1.10 ", \
#     '3Dvlt_K0cFind2_Plot_Cc1.00_try3.dat' using 2:3 with points title " '2try3' 1.00 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.90_try3.dat' using 2:3 with points title " '2try3' 0.90 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.80_try3.dat' using 2:3 with points title " '2try3' 0.80 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.10.dat' using 2:3 with points title " '2' 0.10 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.20.dat' using 2:3 with points title " '2' 0.20 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.30.dat' using 2:3 with points title " '2' 0.30 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.40.dat' using 2:3 with points title " '2' 0.40 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.50.dat' using 2:3 with points title " '2' 0.50 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.60.dat' using 2:3 with points title " '2' 0.60 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.70.dat' using 2:3 with points title " '2' 0.70 ", \
#     '3Dvlt_K0cFind2_Plot_Cc1.20.dat' using 2:3 with points title " '2' 1.20 ", \
#     '3Dvlt_K0cFind2_Plot_Cc1.10.dat' using 2:3 with points title " '2' 1.10 ", \
#     '3Dvlt_K0cFind2_Plot_Cc1.00.dat' using 2:3 with points title " '2' 1.00 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.90.dat' using 2:3 with points title " '2' 0.90 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.80.dat' using 2:3 with points title " '2' 0.80 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.10_try.dat' using 2:3 with points title " '2try' 0.10 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.20_try.dat' using 2:3 with points title " '2try' 0.20 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.30_try.dat' using 2:3 with points title " '2try' 0.30 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.40_try.dat' using 2:3 with points title " '2try' 0.40 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.50_try.dat' using 2:3 with points title " '2try' 0.50 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.60_try.dat' using 2:3 with points title " '2try' 0.60 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.70_try.dat' using 2:3 with points title " '2try' 0.70 ", \
#     '3Dvlt_K0cFind2_Plot_Cc1.20_try.dat' using 2:3 with points title " '2try' 1.20 ", \
#     '3Dvlt_K0cFind2_Plot_Cc1.10_try.dat' using 2:3 with points title " '2try' 1.10 ", \
#     '3Dvlt_K0cFind2_Plot_Cc1.00_try.dat' using 2:3 with points title " '2try' 1.00 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.90_try.dat' using 2:3 with points title " '2try' 0.90 ", \
#     '3Dvlt_K0cFind2_Plot_Cc0.80_try.dat' using 2:3 with points title " '2try' 0.80 ", \
#     '3Dvlt_K0cFind_Plot_Cc0.10_last.dat' using 2:3 with points title " 0.10 ", \
#     '3Dvlt_K0cFind_Plot_Cc0.20_last.dat' using 2:3 with points title " 0.20 ", \
#     '3Dvlt_K0cFind_Plot_Cc0.30_last.dat' using 2:3 with points title " 0.30 ", \
#     '3Dvlt_K0cFind_Plot_Cc0.40_last.dat' using 2:3 with points title " 0.40 ", \
#     '3Dvlt_K0cFind_Plot_Cc0.50_last.dat' using 2:3 with points title " 0.50 ", \
#     '3Dvlt_K0cFind_Plot_Cc0.60_last.dat' using 2:3 with points title " 0.60 ", \
#     '3Dvlt_K0cFind_Plot_Cc0.70_last.dat' using 2:3 with points title " 0.70 ", \
#     '3Dvlt_K0cFind_Plot_Cc0.80_last.dat' using 2:3 with points title " 0.80 ", \
#     '3Dvlt_K0cFind_Plot_Cc0.90_last.dat' using 2:3 with points title " 0.90 ", \
#     '3Dvlt_K0cFind_Plot_Cc1.00_last.dat' using 2:3 with points title " 1.00 ", \
#     '3Dvlt_K0cFind_Plot_Cc1.10_last.dat' using 2:3 with points title " 1.10 ", \
#     '3Dvlt_K0cFind_Plot_Cc1.20_last.dat' using 2:3 with points title " 1.20 ", \





#     "<echo 0.38750818971 0.062421005458"
#     'FixedPoints.dat' using 2:3 with points lc rgb "black" title " Fixed pts "



#ls 2 lc rgb "red"     title " Adjusted theory, {/Times-Italic c}@_{/Times-Italic P}^{adj} ", \
#ls 2 lc rgb "blue"    notitle, \
#ls 2 lc rgb "magenta" notitle, \
#ls 2 lc rgb "cyan"    notitle, \
#ls 2 lc rgb "yellow"  notitle, \
#ls 2 lc rgb "black"   notitle

#'3Dvlt_K0cFind_Plot_Cc0.10_K015.dat' using 2:3 with points title " 0.10 ", \
#     '3Dvlt_K0cFind_Plot_Cc0.20_K015.dat' using 2:3 with points title " 0.20 ", \
#     '3Dvlt_K0cFind_Plot_Cc0.30_K015.dat' using 2:3 with points title " 0.30 ", \
#     '3Dvlt_K0cFind_Plot_Cc0.40_K015.dat' using 2:3 with points title " 0.40 ", \
#     '3Dvlt_K0cFind_Plot_Cc0.50_K015.dat' using 2:3 with points title " 0.50 ", \
#     '3Dvlt_K0cFind_Plot_Cc0.60_K015.dat' using 2:3 with points title " 0.60 ", \
#     '3Dvlt_K0cFind_Plot_Cc0.70_K015.dat' using 2:3 with points title " 0.70 ", \
#     '3Dvlt_K0cFind_Plot_Cc0.80_K015.dat' using 2:3 with points title " 0.80 ", \
#     '3Dvlt_K0cFind_Plot_Cc0.90_K015.dat' using 2:3 with points title " 0.90 ", \
#     '3Dvlt_K0cFind_Plot_Cc1.00_K015.dat' using 2:3 with points title " 1.00 ", \
#     '3Dvlt_K0cFind_Plot_Cc1.10_K015.dat' using 2:3 with points title " 1.10 ", \
#     '3Dvlt_K0cFind_Plot_Cc1.20_K015.dat' using 2:3 with points title " 1.20 ", \
