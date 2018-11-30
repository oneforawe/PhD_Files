#filename: thesis_plot_3Dvlt_flow_itlog_UokT.gnu
reset
set terminal postscript color eps enhanced
set output 'thesis_plot_3Dvlt_flow_itlog_UokT.ps'
set size 0.85,0.85
##set size 1.00,0.85
#set size 0.85,1.20

#set samples 1000000

set title "^4He Vortex Energy {/Times-Italic U/kT} vs Length Scale Parameter {/Times-Italic l}"

set ylabel "{/Times-Italic U/kT} (unitless)"
set xlabel "{/Times-Italic l} (unitless)"



PISQ = 9.86960440108935861883
PICU = 31.00627668029981620634
A3 = 4.0*PICU/3.0

KstarD = 0.607927101854026652
ystarD = 0.039788735772973836

KstarA = 4.146766416946759293
ystarA = 0.0058331356032128735

KstarB = 0.387508189712343409
ystarB = 0.0624210054576019

#peek
#set xrange [-0.1:12.0]
#set yrange [-0.01:0.6]


#set log y
#set xrange [1.0e-10:1.0e5]
#set xrange [-5.0:250.0]
#set xrange [-5.0:5000.0]
#set yrange [1.0e-323:1.0e10]
#set xrange [-5.0:35.0]


#set xrange [-0.1:10.0]
#set xrange [-0.1:8.0]
###############################set xrange [-0.1:5.0]
###############################set yrange [-0.01:0.26]

#set xrange [KstarD-0.25:KstarD+0.25]
#set yrange [ystarD-0.025:ystarD+0.025]

#set xrange [6:8]
#set yrange [0:0.01]
#set xrange [-0.1:50.0]
#set yrange [-0.4:0.4]
#set xrange [-0.006:0.006]
#set yrange [-0.002:0.002]

#set xrange [-0.1:50]
#set yrange [-0.1:0.4]
#set xrange [-0.1:35]
#set yrange [-0.01:0.4]

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
#set xrange [KstarD-3e-4:KstarD+3e-4]
#set yrange [ystarD-5e-5:ystarD+5e-5]

#set xrange [KstarD-3e-5:KstarD+3e-5]
#set yrange [ystarD-5e-6:ystarD+5e-6]

###set xrange [KstarD-3e-6:KstarD+3e-6]
###set yrange [ystarD-5e-7:ystarD+5e-7]
####set xrange [KstarD-1.1e-6:KstarD+1.1e-6]
####set yrange [ystarD-4e-7:ystarD+4e-7]

#set xrange [KstarD-3e-7:KstarD+3e-7]
#set yrange [ystarD-5e-8:ystarD+5e-8]

#set xrange [KstarD-3e-14:KstarD+3e-14]
#set yrange [ystarD-5e-15:ystarD+5e-15]
#set xrange [KstarD-3e-15:KstarD+3e-15]
#set yrange [ystarD-5e-16:ystarD+5e-16]
#set xrange [KstarD-3e-16:KstarD+3e-16]
#set yrange [ystarD-5e-17:ystarD+5e-17]

##set xrange [3.4:5.2]
##set yrange [0.003:0.012]

set grid
set key inside top right box width -5 height 0.5 spacing 1.2 title " ^{}Bare ({/Times-Italic l} = 0) values "
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


plot '3Dvlt_flow_Op11133_K0.212_y0.100_step1000000.dat' using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.141_y0.216_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.156_y0.184_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.174_y0.151_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.198_y0.117_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.247_y0.068_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.282_y0.047_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.296_y0.040_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.312_y0.034_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.371_y0.018_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.494_y0.005_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.741_y0.000_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.262_y0.356_step1000000.dat' using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.288_y0.321_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.320_y0.283_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.360_y0.242_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.411_y0.197_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.479_y0.151_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.548_y0.115_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.575_y0.103_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.606_y0.092_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.719_y0.058_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K0.959_y0.023_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K1.438_y0.003_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle, \
     '3Dvlt_flow_Op11133_K2.877_y0.000_step1000000.dat'      using 4:(6.0*$4-log($3)) with lines ls 2 lc rgb "black"   notitle

#1.0/(A3*x) notitle ls 0, \
#     '3Dvlt_flow_plot_Op13264_K0.388_y0.062_step1000000.dat' using 2:3 with lines ls 1 lc rgb "red"     notitle, \
#     '3Dvlt_flow_plot_Op13254_K0.388_y0.062_step1000000.dat' using 2:3 with lines ls 1 lc rgb "blue"    notitle, \
#     '3Dvlt_flow_plot_Op13363_K0.388_y0.062_step1000000.dat' using 2:3 with lines ls 1 lc rgb "green"   notitle, \
#     '3Dvlt_flow_plot_Op13353_K0.388_y0.062_step1000000.dat' using 2:3 with lines ls 1 lc rgb "magenta" notitle, \
#     exp(-PISQ*0.40*x) lt 5 lw 1.5 lc rgb "orange" title "C = 0.4", \
#     exp(-PISQ*1.10*x) lt 6 lw 1.5 lc rgb "purple" title "C = 1.1", \
#     "<echo 0.607927101854026652 0.039788735772973836"  with points ls 7 lc rgb "#888888" notitle, \
#     "<echo 0.0                  0.0                 "  with points ls 7 lc rgb "#888888" notitle, \
#     "<echo 0.387508189712343409 0.0624210054576019"    with points ls 7 lc rgb "black" notitle, \
#     "<echo 4.146766416946759293 0.0058331356032128735" with points ls 7 lc rgb "#888888" notitle
