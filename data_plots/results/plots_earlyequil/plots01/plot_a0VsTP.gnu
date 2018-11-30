#filename: plot_a0VsTP.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_a0VsTP.ps'
set size 1.3,1

set title "^4He Minimal Vortex Diameter {/Times-Italic a}_0 vs Pressure {/Times-Italic P} \n at Different Temperatures {/Times-Italic T} "

set ylabel "{/Times-Italic a}_0 (angstrom)"
set xlabel "{/Times-Italic P} (bar)"
#set logscale x
set xrange [0:24]
set grid
set key outside right box width -10 spacing 1.2 title " Data from {/CM-Typewriter RhosAmps.ods} "
set rmargin 17

plot 'a0VsTP.dat' using 1:2  with linespoints title " {/Times-Italic T} = 0.0 K ", \
     'a0VsTP.dat' using 1:3  with linespoints title " {/Times-Italic T} = 0.1 K ", \
     'a0VsTP.dat' using 1:4  with linespoints title " {/Times-Italic T} = 0.2 K ", \
     'a0VsTP.dat' using 1:5  with linespoints title " {/Times-Italic T} = 0.3 K ", \
     'a0VsTP.dat' using 1:6  with linespoints title " {/Times-Italic T} = 0.4 K ", \
     'a0VsTP.dat' using 1:7  with linespoints title " {/Times-Italic T} = 0.5 K ", \
     'a0VsTP.dat' using 1:8  with linespoints title " {/Times-Italic T} = 0.6 K ", \
     'a0VsTP.dat' using 1:9  with linespoints title " {/Times-Italic T} = 0.7 K ", \
     'a0VsTP.dat' using 1:10 with linespoints title " {/Times-Italic T} = 0.8 K ", \
     'a0VsTP.dat' using 1:11 with linespoints title " {/Times-Italic T} = 0.9 K ", \
     'a0VsTP.dat' using 1:12 with linespoints title " {/Times-Italic T} = 1.0 K ", \
     'a0VsTP.dat' using 1:13 with linespoints title " {/Times-Italic T} = 1.1 K ", \
     'a0VsTP.dat' using 1:14 with linespoints title " {/Times-Italic T} = 1.2 K ", \
     'a0VsTP.dat' using 1:15 with linespoints title " {/Times-Italic T} = 1.3 K ", \
     'a0VsTP.dat' using 1:16 with linespoints title " {/Times-Italic T} = 1.4 K ", \
     'a0VsTP.dat' using 1:17 with linespoints title " {/Times-Italic T} = 1.5 K ", \
     'a0VsTP.dat' using 1:18 with linespoints title " {/Times-Italic T} = 1.6 K ", \
     'a0VsTP.dat' using 1:19 with linespoints title " {/Times-Italic T} = 1.7 K ", \
     'a0VsTP.dat' using 1:20 with linespoints title " {/Times-Italic T} = 1.8 K ", \
     'a0VsTP.dat' using 1:21 with linespoints title " {/Times-Italic T} = 1.9 K ", \
     'a0VsTP.dat' using 1:22 with linespoints title " {/Times-Italic T} = 2.0 K ", \
     'a0VsTP.dat' using 1:23 with linespoints title " {/Times-Italic T} = 2.1 K ", \
     'a0VsTP.dat' using 1:24 with linespoints title " {/Times-Italic T} = 2.2 K ", \
     'a0VsTP.dat' using 1:25 with linespoints title " {/Times-Italic T} = 2.3 K ", \
     'a0VsTP.dat' using 1:26 with linespoints title " {/Times-Italic T} = 2.4 K ", \
     'a0VsTP.dat' using 1:27 with linespoints title " {/Times-Italic T} = 2.5 K "
