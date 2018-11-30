#filename: thesis_plot_3Dvlt_K0cFind3_Plot_Ky_full.gnu
reset
set terminal postscript color eps enhanced
set output 'thesis_plot_3Dvlt_K0cFind3_Plot_Ky_full.ps'
set size 0.85,0.85
##set size 1.00,0.85
#set size 0.85,1.20

#set samples 1000000

set title "^4He Fugacity {/Times-Italic y} vs Coupling Parameter {/Times-Italic K}"

set ylabel "{/Times-Italic y} (unitless)"
set xlabel "{/Times-Italic K} (unitless)"

#set log y
#set xrange [1.0e-10:1.0e5]
#set xrange [-5.0:250.0]
#set xrange [-5.0:5000.0]
#set yrange [1.0e-323:1.0e10]
#set xrange [-5.0:35.0]



#set xrange [0.387508189712343409-0.25:0.387508189712343409+0.25]
#set yrange [0.0624210054576019-0.025:0.0624210054576019+0.025]

#set xrange [6:8]
#set yrange [0:0.01]
set xrange [-0.1:35.0]
set yrange [-0.4:0.4]
##set yrange [-0.025:0.025]
#set xrange [-0.006:0.006]
#set yrange [-0.002:0.002]

#set xrange [-0.1:50]
#set yrange [-0.1:0.4]
#set xrange [-0.1:35]
#set yrange [-0.01:0.4]
#set xrange [-0.1:5]
#set yrange [-0.01:0.1]

#set xrange [30.5:31.5]
#set yrange [-0.01:0.01]

#set xrange [-0.1:32]
#set yrange [-0.1:0.3]
##set xrange [-0.1:9]
##set yrange [-0.03:0.27]
##############set xrange [0.2:2]
##############set yrange [0.02:0.22]
#set xrange [-0.1:9]
#set yrange [-0.03:0.1]
#set xrange [-0.1:10.5]
#set yrange [-0.01:0.1]

#set xrange [-0.1:5.5]
#set yrange [-0.01:0.3]
#set xrange [-0.1:1]
#set yrange [-0.01:0.1]

##set xrange [0.25:0.5]
##set yrange [0.04:0.08]
###set xrange [0.2:1.5]
###set yrange [0.01:0.15]
#set xrange [0.387508189712343409-2e-4:0.387508189712343409+2e-4]
#set yrange [0.0624210054576019-5e-5:0.0624210054576019+5e-5]

#set xrange [0.387508189712343409-3e-5:0.387508189712343409+3e-5]
#set yrange [0.0624210054576019-5e-6:0.0624210054576019+5e-6]

###set xrange [0.387508189712343409-3e-6:0.387508189712343409+3e-6]
###set yrange [0.0624210054576019-5e-7:0.0624210054576019+5e-7]
####set xrange [0.387508189712343409-1.1e-6:0.387508189712343409+1.1e-6]
####set yrange [0.0624210054576019-4e-7:0.0624210054576019+4e-7]

#set xrange [0.387508189712343409-3e-7:0.387508189712343409+3e-7]
#set yrange [0.0624210054576019-5e-8:0.0624210054576019+5e-8]

#set xrange [0.387508189712343409-3e-14:0.387508189712343409+3e-14]
#set yrange [0.0624210054576019-5e-15:0.0624210054576019+5e-15]
#set xrange [0.387508189712343409-3e-15:0.387508189712343409+3e-15]
#set yrange [0.0624210054576019-5e-16:0.0624210054576019+5e-16]
#set xrange [0.387508189712343409-3e-16:0.387508189712343409+3e-16]
#set yrange [0.0624210054576019-5e-17:0.0624210054576019+5e-17]

##set xrange [3.4:5.2]
##set yrange [0.003:0.012]

set grid
#####set key at graph 0.85,-0.21 box width -5 height 0.5 spacing 1.6 title " Finding {/Times-Italic K}@^{{} {/Times-Italic l}0^{^{}}}_{c  }, given a {/Times-Italic C}_c^ "
##set key at graph 0.875,-0.21 box width -5 height 0.5 spacing 1.6 title " Finding {/Times-Italic K}@^{{} {/Times-Italic l}0^{^{}}}_{c  }, given a {/Times-Italic C}_c^ "
#set key outside below center box width -5 height 0.5 spacing 1.2 title " Finding {/Times-Italic K}_{0c}, given a {/Times-Italic C}@_c^ "
#####set bmargin 9.7

set rmargin 5
##set key outside right
#set bmargin 20
#set key outside below
set key off

#set pointsize 0.4
set pointsize 0.5
#set pointsize 0.125

PISQ = 9.86960440108935861883
PICU = 31.00627668029981620634
A3 = 4.0*PICU/3.0


plot '3Dvlt_K0cFind3_Plot_Op3341_K32.600_y0.390_step1000000.dat'  using 2:3 with lines ls 2 lc rgb "black" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3352_K32.600_y0.390_step1000000.dat'  using 2:3 with lines ls 2 lc rgb "black" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3341_K32.600_y0.290_step1000000.dat'  using 2:3 with lines ls 2 lc rgb "black" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3352_K32.600_y0.290_step1000000.dat'  using 2:3 with lines ls 2 lc rgb "black" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3341_K32.600_y0.190_step1000000.dat'  using 2:3 with lines ls 2 lc rgb "black" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3352_K32.600_y0.190_step1000000.dat'  using 2:3 with lines ls 2 lc rgb "black" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3341_K32.600_y0.090_step1000000.dat'  using 2:3 with lines ls 2 lc rgb "black" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3352_K32.600_y0.090_step1000000.dat'  using 2:3 with lines ls 2 lc rgb "black" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3341_K32.600_y-0.010_step1000.dat'    using 2:3 with lines ls 2 lc rgb "orange" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3352_K32.600_y-0.010_step1000000.dat' using 2:3 with lines ls 2 lc rgb "orange" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3341_K32.600_y-0.110_step1000.dat'    using 2:3 with lines ls 2 lc rgb "orange" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3352_K32.600_y-0.110_step1000000.dat' using 2:3 with lines ls 2 lc rgb "orange" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3341_K32.600_y-0.210_step1000.dat'    using 2:3 with lines ls 2 lc rgb "orange" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3352_K32.600_y-0.210_step1000000.dat' using 2:3 with lines ls 2 lc rgb "orange" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3341_K32.600_y-0.310_step1000.dat'    using 2:3 with lines ls 2 lc rgb "orange" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3352_K32.600_y-0.310_step1000000.dat' using 2:3 with lines ls 2 lc rgb "orange" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3341_K32.600_y-0.410_step1000.dat'    using 2:3 with lines ls 2 lc rgb "orange" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3352_K32.600_y-0.410_step1000000.dat' using 2:3 with lines ls 2 lc rgb "orange" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3341_K32.600_y-0.510_step1000.dat'    using 2:3 with lines ls 2 lc rgb "orange" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3352_K32.600_y-0.510_step1000000.dat' using 2:3 with lines ls 2 lc rgb "orange" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3352_K32.600_y-0.610_step1000000.dat' using 2:3 with lines ls 2 lc rgb "orange" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3363_K0.100_y-0.001_step4000000.dat'  using 2:3 with lines ls 2 lc rgb "cyan" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3341_K0.388_y-0.075_step1000000.dat'  using 2:3 with lines ls 2 lc rgb "cyan" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3352_K0.388_y-0.075_step1000000.dat'  using 2:3 with lines ls 2 lc rgb "cyan" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3341_K4.145_y0.157_step100000.dat'    using 2:3 with lines ls 2 lc rgb "black" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3352_K4.145_y0.157_step100000.dat'    using 2:3 with lines ls 2 lc rgb "black" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3341_K4.145_y0.057_step1000000.dat'   using 2:3 with lines ls 2 lc rgb "black" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op3352_K4.145_y0.057_step100000.dat'    using 2:3 with lines ls 2 lc rgb "black" title "(ds limtd)", \
     '3Dvlt_K0cFind3_Plot_Op2134_K0.388_y0.062_step11000000.dat' using 2:3 with lines ls 1 lc rgb "red"     title "(dl const)", \
     '3Dvlt_K0cFind3_Plot_Op2122_K0.388_y0.062_step4000000.dat'  using 2:3 with lines ls 1 lc rgb "blue"    title "(dl const)", \
     '3Dvlt_K0cFind3_Plot_Op2231_K0.388_y0.062_step2000000.dat'  using 2:3 with lines ls 1 lc rgb "green"   title "(dl const)", \
     '3Dvlt_K0cFind3_Plot_Op2221_K0.388_y0.062_step1500000.dat'  using 2:3 with lines ls 1 lc rgb "magenta" title "(dl const)", \
     "<echo 0.387508189712343409409 0.0624210054576019" with points ls 7 lc rgb "black" title "Fixed pts", \
     "<echo 4.146766416946759293 0.0058331356032128735" with points ls 7 lc rgb "black" notitle


#     '3Dvlt_K0cFind3_Plot_Op3341_K1.000_y0.000_step1500000.dat'  using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3352_K1.000_y0.000_step1000000.dat'  using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \

#     '3Dvlt_K0cFind3_Plot_Op3352_K32.600_y-0.710_step1000000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(ds limtd)", \


#'3Dvlt_K0cFind3_Plot_Op3363_K0.005_y-0.001_step4000000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(ds limtd)", \
#     '3Dvlt_K0cFind3_Plot_Op3363_K0.010_y-0.001_step4000000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(ds limtd)", \
##     '3Dvlt_K0cFind3_Plot_Op3363_K0.040_y-0.001_step4000000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(ds limtd)", \
##     '3Dvlt_K0cFind3_Plot_Op3363_K0.050_y-0.001_step4000000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(ds limtd)", \
##     '3Dvlt_K0cFind3_Plot_Op3363_K0.070_y-0.001_step4000000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(ds limtd)", \


#     '3Dvlt_K0cFind3_Plot_Op3341_K0.38755_y0.06242_step100000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3352_K0.38755_y0.06242_step100000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \

#     '3Dvlt_K0cFind3_Plot_Op3341_K25.000_y0.000_step1000000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3352_K25.000_y0.000_step1000000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3352_K32.000_y0.000_step1000000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3352_K40.000_y0.000_step100000.dat'  using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \

#     '3Dvlt_K0cFind3_Plot_Op2122_K0.388_y0.062_step1000000.dat'  using 2:3 with lines ls 1 lc rgb "blue"    title "(dl const)", \




#     '3Dvlt_K0cFind3_Plot_Op3341_K0.38750_y0.06240_step100000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3352_K0.38750_y0.06240_step100000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3341_K0.38745_y0.06242_step150000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3352_K0.38745_y0.06242_step100000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3341_K0.38758_y0.06242_step150000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3352_K0.38758_y0.06242_step100000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3341_K0.38750_y0.06244_step100000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3352_K0.38750_y0.06244_step100000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3341_K0.38738_y0.06244_step100000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3352_K0.38738_y0.06244_step100000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3341_K0.38750_y0.06238_step100000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3352_K0.38750_y0.06238_step100000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3341_K0.38750_y0.06246_step100000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3352_K0.38750_y0.06246_step100000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3341_K0.38765_y0.06240_step100000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3352_K0.38765_y0.06240_step100000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \



#     '3Dvlt_K0cFind3_Plot_Op3341_K0.400_y0.025_step150000.dat'    using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3352_K0.400_y0.025_step1200000.dat'  using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3341_K0.400_y0.090_step500000.dat'   using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3352_K0.400_y0.090_step1200000.dat'  using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3341_K40.000_y0.000_step200000.dat'  using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3352_K40.000_y0.000_step1000000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3341_K100.000_y0.000_step500000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3352_K100.000_y0.000_step100000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3352_K1.000_y1.000_step50000000.dat' using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3341_K1.000_y0.000_step1500000.dat'  using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3352_K1.000_y0.000_step1000000.dat'  using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3341_K0.000_y0.000_step1500000.dat'  using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3352_K1.000_y3.000_step500000.dat'   using 2:3 with lines ls 2 lc rgb "black" title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3341_K0.000_y0.000_step1000000.dat'  using 2:3 with lines ls 2 lc rgb "black"   title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3341_K0.150_y0.060_step300000.dat'   using 2:3 with lines ls 2 lc rgb "black"   title "(dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op3352_K0.150_y0.060_step1000000.dat'  using 2:3 with lines ls 2 lc rgb "black"   title "(dl const)", \






#     '3Dvlt_K0cFind3_Plot_Op3374_K0.250_y0.250_step11000000.dat' using 2:3 with points title "(ds limtd)", \
#     '3Dvlt_K0cFind3_Plot_Op3363_K0.250_y0.013_step4000000.dat' using 2:3 with points title "(ds limtd)", \
#     '3Dvlt_K0cFind3_Plot_Op3363_K3.000_y0.050_step4000000.dat' using 2:3 with points title "(ds limtd)", \
#     '3Dvlt_K0cFind3_Plot_Op3374_K3.000_y0.050_step4000000.dat' using 2:3 with points title "(ds limtd)", \
     
#     exp(-PISQ*0.70*x) lt 2 notitle, \
#     exp(-PISQ*0.71*x) lt 2 notitle, \
#     exp(-PISQ*0.72*x) lt 2 notitle, \
#     exp(-PISQ*0.73*x) lt 2 notitle, \
#     exp(-PISQ*0.74*x) lt 2 notitle, \
#     exp(-PISQ*0.75*x) lt 2 notitle, \



#'3Dvlt_K0cFind3_Plot_2114Op_Cc0.10.dat' using 2:3 with points title "Cc = 0.10 (ds limtd)", \
#     '3Dvlt_K0cFind3_Plot_2114Op_Cc0.20.dat' using 2:3 with points title "Cc = 0.20 (ds limtd)", \
#     '3Dvlt_K0cFind3_Plot_2114Op_Cc0.30.dat' using 2:3 with points title "Cc = 0.30 (ds limtd)", \
#     '3Dvlt_K0cFind3_Plot_2114Op_Cc0.40.dat' using 2:3 with points title "Cc = 0.40 (ds limtd)", \
#     '3Dvlt_K0cFind3_Plot_2114Op_Cc0.50.dat' using 2:3 with points title "Cc = 0.50 (ds limtd)", \
#     '3Dvlt_K0cFind3_Plot_2114Op_Cc0.60.dat' using 2:3 with points title "Cc = 0.60 (ds limtd)", \
#     '3Dvlt_K0cFind3_Plot_2114Op_Cc0.70.dat' using 2:3 with points title "Cc = 0.70 (ds limtd)", \
#     '3Dvlt_K0cFind3_Plot_2114Op_Cc1.20.dat' using 2:3 with points title "Cc = 1.20 (ds limtd)", \
#     '3Dvlt_K0cFind3_Plot_2114Op_Cc1.10.dat' using 2:3 with points title "Cc = 1.10 (ds limtd)", \
#     '3Dvlt_K0cFind3_Plot_2114Op_Cc1.00.dat' using 2:3 with points title "Cc = 1.00 (ds limtd)", \
#     '3Dvlt_K0cFind3_Plot_2114Op_Cc0.90.dat' using 2:3 with points title "Cc = 0.90 (ds limtd)", \
#     '3Dvlt_K0cFind3_Plot_2114Op_Cc0.80.dat' using 2:3 with points title "Cc = 0.80 (ds limtd)", \
#     '3Dvlt_K0cFind3_Plot_2112Op_Cc0.10.dat' using 2:3 with points title "Cc = 0.10 (dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op2112_Cc0.40_uncert1e-08.dat' using 2:3 with points title "Cc = 0.40 (dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op2112_Cc0.50_uncert1e-08.dat' using 2:3 with points title "Cc = 0.50 (dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op2112_Cc0.60_uncert1e-08.dat' using 2:3 with points title "Cc = 0.60 (dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op2112_Cc0.70_uncert1e-08.dat' using 2:3 with points title "Cc = 0.70 (dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op2112_Cc1.20_uncert1e-08.dat' using 2:3 with points title "Cc = 1.20 (dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op2112_Cc1.10_uncert1e-08.dat' using 2:3 with points title "Cc = 1.10 (dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op2112_Cc1.00_uncert1e-08.dat' using 2:3 with points title "Cc = 1.00 (dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op2112_Cc0.90_uncert1e-08.dat' using 2:3 with points title "Cc = 0.90 (dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op2112_Cc0.80_uncert1e-08.dat' using 2:3 with points title "Cc = 0.80 (dl const)", \
#     exp(-PISQ*0.10*x) lt 2 lc rgb "yellow" notitle, \
#     exp(-PISQ*0.20*x) lt 2 lc rgb "yellow" notitle, \
#     exp(-PISQ*0.30*x) lt 2 lc rgb "yellow" notitle, \
#     exp(-PISQ*0.40*x) lt 2 lc rgb "blue" notitle, \
#     exp(-PISQ*0.50*x) lt 2 notitle, \
#     exp(-PISQ*0.60*x) lt 2 notitle, \
#     exp(-PISQ*0.70*x) lt 2 notitle, \
#     exp(-PISQ*0.80*x) lt 2 notitle, \
#     exp(-PISQ*0.90*x) lt 2 notitle, \
#     exp(-PISQ*1.00*x) lt 2 notitle, \
#     exp(-PISQ*1.10*x) lt 2 lc rgb "red"  notitle, \
#     exp(-PISQ*1.20*x) lt 2 lc rgb "yellow" notitle, \
#     1.0/(A3*x) notitle ls 0, \







#    '3Dvlt_K0cFind3_Plot_Op2132_K0.388_y0.062_step2000000.dat' using 2:3 with points title "Cc = 0.80 (dl const)", \
#    '3Dvlt_K0cFind3_Plot_Op2134_K0.388_y0.062_step2000000.dat' using 2:3 with points title "Cc = 0.80 (dl const)", \

#     "<echo 0.387508189712343409 0.0624210054576019" with points lc rgb "black" title " Fixed pts ", \
#     "<echo 0.38750818971 0.062421005458"
#     'FixedPoints.dat' using 2:3 with points lc rgb "black" title " Fixed pts "

#ls 2 lc rgb "red"     title " Adjusted theory, {/Times-Italic c}@_{/Times-Italic P}^{adj} ", \
#ls 2 lc rgb "blue"    notitle, \
#ls 2 lc rgb "magenta" notitle, \
#ls 2 lc rgb "cyan"    notitle, \
#ls 2 lc rgb "yellow"  notitle, \
#ls 2 lc rgb "black"   notitle

#    '3Dvlt_K0cFind3_Plot_3341Op.dat' using 2:3 with points title "", \
#     '3Dvlt_K0cFind3_Plot_3341Op_K7.500_y0.002_step_100000.dat' using 2:3 with points title "", \

#     '3Dvlt_K0cFind3_Plot_Op2112_Cc0.10.dat' using 2:3 with points title "Cc = 0.10 (dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op2112_Cc0.20.dat' using 2:3 with points title "Cc = 0.20 (dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op2112_Cc0.30.dat' using 2:3 with points title "Cc = 0.30 (dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op2112_Cc0.40.dat' using 2:3 with points title "Cc = 0.40 (dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op2112_Cc0.50.dat' using 2:3 with points title "Cc = 0.50 (dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op2112_Cc0.60.dat' using 2:3 with points title "Cc = 0.60 (dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op2112_Cc0.70.dat' using 2:3 with points title "Cc = 0.70 (dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op2112_Cc1.20.dat' using 2:3 with points title "Cc = 1.20 (dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op2112_Cc1.10.dat' using 2:3 with points title "Cc = 1.10 (dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op2112_Cc1.00.dat' using 2:3 with points title "Cc = 1.00 (dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op2112_Cc0.90.dat' using 2:3 with points title "Cc = 0.90 (dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op2112_Cc0.80.dat' using 2:3 with points title "Cc = 0.80 (dl const)", \


#plot '3Dvlt_K0cFind3_Plot_Op2134_K0.388_y0.062_step6000000.dat' using 1:5 with points title "Cc = 0.80 (dl const)", \
#     '3Dvlt_K0cFind3_Plot_Op2134_K0.388_y0.062_step2000000.dat' using 1:5 with points title "Cc = 0.80 (dl const)"
