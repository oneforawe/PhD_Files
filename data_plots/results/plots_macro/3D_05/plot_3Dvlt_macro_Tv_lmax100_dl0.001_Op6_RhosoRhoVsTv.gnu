#filename: plot_3Dvlt_macro_Tv_lmax100_dl0.001_Op6_RhosoRhoVsTv.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_3Dvlt_macro_Tv_lmax100_dl0.001_Op6_RhosoRhoVsTv.ps'
set size 1.3,1.1

set title "^4He Superfluid Fraction {/Symbol-Oblique r}_s / {/Symbol-Oblique r} \n vs Exponentiated Temperature Variable {/Symbol-Oblique t}^{/Symbol-Oblique n}"

set ylabel "{/Symbol-Oblique r}_s / {/Symbol-Oblique r} (unitless)"
set xlabel "{/Symbol-Oblique t} = [{/Times-Italic 1 - T / T}_{/Symbol-Oblique l}{/Times-Italic (P)}] (unitless)"
#set xrange [0:1e-5]
set xrange [0:1e-8]
#set xrange [0:4e-8]
#set yrange [0:3e-5]
#set yrange [0:0.9]
set grid
set key outside right box width -5 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter 3Dvlt\\_macro.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n}"
set rmargin 30

plot '3Dvlt_macro_Tv_Cc0.55_lmax100_dl0.001_Op6.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.55 ", \
     '3Dvlt_macro_Tv_Cc0.50_lmax100_dl0.001_Op6.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.50 ", \
     '3Dvlt_macro_Tv_Cc0.45_lmax100_dl0.001_Op6.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.45 ", \
     '3Dvlt_macro_Tv_Cc0.40_lmax100_dl0.001_Op6.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.40 ", \
     '3Dvlt_macro_Tv_Cc0.35_lmax100_dl0.001_Op6.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.35 ", \
     '3Dvlt_macro_Tv_Cc0.30_lmax100_dl0.001_Op6.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.30 ", \
     '3Dvlt_macro_Tv_Cc0.25_lmax100_dl0.001_Op6.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.25 ", \
     '3Dvlt_macro_Tv_Cc0.20_lmax100_dl0.001_Op6.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.20 ", \
     '3Dvlt_macro_Tv_Cc0.15_lmax100_dl0.001_Op6.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.15 ", \
     '3Dvlt_macro_Tv_Cc0.10_lmax100_dl0.001_Op6.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.10 ", \
     2.64727403089064*x**0.67169157167077  title " 2.64727*{/Symbol-Oblique t}^{0.671692} ", \
     2.48138441221253*x**0.671688975904562 title " 2.48138*{/Symbol-Oblique t}^{0.671689} ", \
     2.28603873457402*x**0.666660017939919 title " 2.28604*{/Symbol-Oblique t}^{0.666660} ", \
     2.52833834313342*x**0.674664792749541 title " 2.52834*{/Symbol-Oblique t}^{0.674665} ", \
     2.60282434534289*x**0.677618255233688 title " 2.60282*{/Symbol-Oblique t}^{0.677618} ", \
     2.35587216841892*x**0.670962750518098 title " 2.35587*{/Symbol-Oblique t}^{0.670963} ", \
     2.3925018679873*x**0.675643069470319  title " 2.39250*{/Symbol-Oblique t}^{0.675643} ", \
     2.11807866887045*x**0.670899109035549 title " 2.11808*{/Symbol-Oblique t}^{0.670899} ", \
     2.04976531298717*x**0.677273734682563 title " 2.04977*{/Symbol-Oblique t}^{0.677274} ", \
     1.6279328232016*x**0.67351783749462   title " 1.62793*{/Symbol-Oblique t}^{0.673518} "
