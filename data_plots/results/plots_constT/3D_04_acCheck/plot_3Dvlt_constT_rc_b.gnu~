#filename: plot_3Dvlt_constT_rc.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_3Dvlt_constT_rc.ps'
set size 0.85,0.85
##set size 1.00,0.85
#set size 0.85,1.20

#set samples 1000000

set title "^4He {/Times-Italic a}_c/{/Times-Italic a} vs Length Scale Parameter {/Times-Italic l}"

set ylabel "{/Times-Italic a}_c/{/Times-Italic a} (unitless)"
set xlabel "{/Times-Italic l} (unitless)"


Kstar = 0.387508189712343409
ystar = 0.0624210054576019

Kstar1 = 4.146766416946759293
ystar1 = 0.0058331356032128735


set xrange [0:110]
set yrange [-0.1:1.1]

#set xrange [Kstar-0.25:Kstar+0.25]
#set yrange [ystar-0.025:ystar+0.025]


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

set pointsize 0.4
#set pointsize 0.5
#set pointsize 0.125

PISQ = 9.86960440108935861883
PICU = 31.00627668029981620634
A3 = 4.0*PICU/3.0


plot '3Dvlt_constT_0.100_Op11_Cc0.40_lmax100_dl1e-05.dat' using 5:13 with lines, \
     '3Dvlt_constT_0.200_Op11_Cc0.40_lmax100_dl1e-05.dat' using 5:13 with lines, \
     '3Dvlt_constT_0.300_Op11_Cc0.40_lmax100_dl1e-05.dat' using 5:13 with lines, \
     '3Dvlt_constT_0.400_Op11_Cc0.40_lmax100_dl1e-05.dat' using 5:13 with lines, \
     '3Dvlt_constT_0.500_Op11_Cc0.40_lmax100_dl1e-05.dat' using 5:13 with lines, \
     '3Dvlt_constT_0.600_Op11_Cc0.40_lmax100_dl1e-05.dat' using 5:13 with lines, \
     '3Dvlt_constT_0.700_Op11_Cc0.40_lmax100_dl1e-05.dat' using 5:13 with lines, \
     '3Dvlt_constT_0.800_Op11_Cc0.40_lmax100_dl1e-05.dat' using 5:13 with lines, \
     '3Dvlt_constT_0.900_Op11_Cc0.40_lmax100_dl1e-05.dat' using 5:13 with lines, \
     '3Dvlt_constT_1.000_Op11_Cc0.40_lmax100_dl1e-05.dat' using 5:13 with lines lc rgb "orange", \
     '3Dvlt_constT_0.100_Op11_Cc1.10_lmax100_dl1e-05.dat' using 5:13 with lines, \
     '3Dvlt_constT_0.200_Op11_Cc1.10_lmax100_dl1e-05.dat' using 5:13 with lines, \
     '3Dvlt_constT_0.300_Op11_Cc1.10_lmax100_dl1e-05.dat' using 5:13 with lines, \
     '3Dvlt_constT_0.400_Op11_Cc1.10_lmax100_dl1e-05.dat' using 5:13 with lines, \
     '3Dvlt_constT_0.500_Op11_Cc1.10_lmax100_dl1e-05.dat' using 5:13 with lines, \
     '3Dvlt_constT_0.600_Op11_Cc1.10_lmax100_dl1e-05.dat' using 5:13 with lines, \
     '3Dvlt_constT_0.700_Op11_Cc1.10_lmax100_dl1e-05.dat' using 5:13 with lines, \
     '3Dvlt_constT_0.800_Op11_Cc1.10_lmax100_dl1e-05.dat' using 5:13 with lines, \
     '3Dvlt_constT_0.900_Op11_Cc1.10_lmax100_dl1e-05.dat' using 5:13 with lines, \
     '3Dvlt_constT_1.000_Op11_Cc1.10_lmax100_dl1e-05.dat' using 5:13 with lines lc rgb "blue"

