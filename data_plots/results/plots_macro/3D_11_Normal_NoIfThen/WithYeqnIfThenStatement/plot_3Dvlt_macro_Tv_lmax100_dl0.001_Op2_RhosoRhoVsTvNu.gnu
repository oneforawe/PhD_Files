#filename: plot_3Dvlt_macro_Tv_lmax100_dl0.001_Op2_RhosoRhoVsTvNu.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_3Dvlt_macro_Tv_lmax100_dl0.001_Op2_RhosoRhoVsTvNu.ps'
set size 1.0,1.1

set title "^4He Superfluid Fraction {/Symbol-Oblique r}_s / {/Symbol-Oblique r} \n vs Exponentiated Temperature Variable {/Symbol-Oblique t}^{/Symbol-Oblique n}"

set ylabel "{/Symbol-Oblique r}_s / {/Symbol-Oblique r} (unitless)"
set xlabel "{/Symbol-Oblique t}^{/Symbol-Oblique n} = [{/Times-Italic 1 - T / T}_{/Symbol-Oblique l}{/Times-Italic (P)}]^{/Symbol-Oblique n} (unitless)"
#set xrange [0:3e-5]
#set yrange [0:1e-4]
set xrange [0:(3.6e-8)**0.671696337558]
set yrange [0:3e-5]
#set xrange [0:4.38e-4]
#set yrange [0:1e-3]
#set xrange [0:4e-5]
#set xrange [0:1e-5]
#set yrange [0:0.1]
set grid
set key outside right box width -16 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter 3Dvlt\\_macro.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.671696337558 "
set rmargin 35

plot '3Dvlt_macro_Tv_Cc1.20_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 1.20 ", \
     '3Dvlt_macro_Tv_Cc1.10_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 1.10 ", \
     '3Dvlt_macro_Tv_Cc1.06_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 1.06 ", \
     '3Dvlt_macro_Tv_Cc1.05_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 1.05 ", \
     '3Dvlt_macro_Tv_Cc1.04_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 1.04 ", \
     '3Dvlt_macro_Tv_Cc1.03_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 1.03 ", \
     '3Dvlt_macro_Tv_Cc1.02_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 1.02 ", \
     '3Dvlt_macro_Tv_Cc1.01_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 1.01 ", \
     '3Dvlt_macro_Tv_Cc1.00_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 1.00 ", \
     '3Dvlt_macro_Tv_Cc0.99_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 0.99 ", \
     '3Dvlt_macro_Tv_Cc0.98_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 0.98 ", \
     '3Dvlt_macro_Tv_Cc0.97_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 0.97 ", \
     '3Dvlt_macro_Tv_Cc0.90_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 0.90 ", \
     '3Dvlt_macro_Tv_Cc0.80_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 0.80 ", \
     '3Dvlt_macro_Tv_Cc0.70_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 0.70 ", \
     '3Dvlt_macro_Tv_Cc0.60_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 0.60 ", \
     '3Dvlt_macro_Tv_Cc0.55_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 0.55 ", \
     '3Dvlt_macro_Tv_Cc0.50_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 0.50 ", \
     '3Dvlt_macro_Tv_Cc0.40_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 0.40 ", \
     '3Dvlt_macro_Tv_Cc0.30_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 0.30 ", \
     '3Dvlt_macro_Tv_Cc0.20_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 0.20 ", \
     4.13069939912516*x title " {/Times-Italic A}{/Symbol \242} = 4.13070 ", \
     3.9233933229044*x title " {/Times-Italic A}{/Symbol \242} = 3.92339 ", \
     3.83483018432706*x title " {/Times-Italic A}{/Symbol \242} = 3.83483 ", \
     3.81215101649793*x title " {/Times-Italic A}{/Symbol \242} = 3.81215 ", \
     3.78925121651428*x title " {/Times-Italic A}{/Symbol \242} = 3.78925 ", \
     3.76612549137784*x title " {/Times-Italic A}{/Symbol \242} = 3.76613 ", \
     3.74277082421197*x title " {/Times-Italic A}{/Symbol \242} = 3.74277 ", \
     3.71918569831501*x title " {/Times-Italic A}{/Symbol \242} = 3.71919 ", \
     3.69536745754452*x title " {/Times-Italic A}{/Symbol \242} = 3.69537 ", \
     3.67130784114546*x title " {/Times-Italic A}{/Symbol \242} = 3.67131 ", \
     3.64701474862901*x title " {/Times-Italic A}{/Symbol \242} = 3.64701 ", \
     3.62247266798771*x title " {/Times-Italic A}{/Symbol \242} = 3.62247 ", \
     3.79854713441012*x title " {/Times-Italic A}{/Symbol \242} = 3.79855 ", \
     3.33084706176453*x title " {/Times-Italic A}{/Symbol \242} = 3.33085 ", \
     3.74327907216237*x title " {/Times-Italic A}{/Symbol \242} = 3.74328 ", \
     2.63831363412165*x title " {/Times-Italic A}{/Symbol \242} = 2.63831 ", \
     3.74805174003843*x title " {/Times-Italic A}{/Symbol \242} = 3.74805 ", \
     3.47257606634826*x title " {/Times-Italic A}{/Symbol \242} = 3.47258 ", \
     3.38447001700223*x title " {/Times-Italic A}{/Symbol \242} = 3.38447 ", \
     2.94851690242369*x title " {/Times-Italic A}{/Symbol \242} = 2.94852 ", \
     1.75806362943721*x title " {/Times-Italic A}{/Symbol \242} = 1.75806 "

#     '3Dvlt_macro_Tv_Cc0.10_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 0.10 ", \
