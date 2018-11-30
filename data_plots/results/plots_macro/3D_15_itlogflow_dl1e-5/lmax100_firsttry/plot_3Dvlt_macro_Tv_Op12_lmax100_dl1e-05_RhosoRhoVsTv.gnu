#filename: plot_3Dvlt_macro_Tv_Op12_lmax100_dl1e-05_RhosoRhoVsTv.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_3Dvlt_macro_Tv_Op12_lmax100_dl1e-05_RhosoRhoVsTv.ps'
set size 1.3,1.1

set title "^4He Superfluid Fraction {/Symbol-Oblique r}_s / {/Symbol-Oblique r} \n vs Exponentiated Temperature Variable {/Symbol-Oblique t}^{/Symbol-Oblique n}"

set ylabel "{/Symbol-Oblique r}_s / {/Symbol-Oblique r} (unitless)"
set xlabel "{/Symbol-Oblique t} = [{/Times-Italic 1 - T / T}_{/Symbol-Oblique l}{/Times-Italic (P)}] (unitless)"
set xrange [0:1e-5]
##set xrange [0:5e-8]
#set xrange [0:4e-8]
#set yrange [0:3e-5]
#set yrange [0:0.9]
set grid
set key outside right box width -13 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter 3Dvlt\\_macro.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n}"
set rmargin 40

plot '3Dvlt_macro_Tv_Op12_Cc0.20_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.20 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.25_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.25 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.30_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.30 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.35_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.35 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.40_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.40 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.45_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.45 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.50_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.50 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.55_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.55 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.60_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.60 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.65_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.65 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.70_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.70 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.75_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.75 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.80_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.80 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.85_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.85 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.90_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.90 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.95_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.95 ", \
     '3Dvlt_macro_Tv_Op12_Cc1.00_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.00 ", \
     '3Dvlt_macro_Tv_Op12_Cc1.05_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 1.05 ", \
     '3Dvlt_macro_Tv_Op12_Cc1.10_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 1.10 ", \
     '3Dvlt_macro_Tv_Op12_Cc1.15_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 1.15 ", \
     '3Dvlt_macro_Tv_Op12_Cc1.20_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 1.20 "
