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
     '3Dvlt_macro_Tv_Cc0.10_lmax100_dl0.001_Op2.dat' using ($2**0.671696337558):6 with points title " {/Times-Italic C}_c = 0.10 ", \
     4.01187402369329*x title " {/Times-Italic A}{/Symbol \242} = 4.01187 ", \
     3.86677170432738*x title " {/Times-Italic A}{/Symbol \242} = 3.86677 ", \
     3.80469569424847*x title " {/Times-Italic A}{/Symbol \242} = 3.80470 ", \
     3.78878357845808*x title " {/Times-Italic A}{/Symbol \242} = 3.78878 ", \
     3.77270486968496*x title " {/Times-Italic A}{/Symbol \242} = 3.77270 ", \
     3.75646856136679*x title " {/Times-Italic A}{/Symbol \242} = 3.75647 ", \
     3.74005717091211*x title " {/Times-Italic A}{/Symbol \242} = 3.74006 ", \
     3.72347695050249*x title " {/Times-Italic A}{/Symbol \242} = 3.72348 ", \
     3.7067158625785*x title " {/Times-Italic A}{/Symbol \242} = 3.70672 ", \
     3.68978659434611*x title " {/Times-Italic A}{/Symbol \242} = 3.68979 ", \
     3.67266976176842*x title " {/Times-Italic A}{/Symbol \242} = 3.67267 ", \
     3.65536737134034*x title " {/Times-Italic A}{/Symbol \242} = 3.65537 ", \
     3.52875068354362*x title " {/Times-Italic A}{/Symbol \242} = 3.52875 ", \
     3.32899509968019*x title " {/Times-Italic A}{/Symbol \242} = 3.32900 ", \
     3.10236821977425*x title " {/Times-Italic A}{/Symbol \242} = 3.10237 ", \
     2.84195365299783*x title " {/Times-Italic A}{/Symbol \242} = 2.84195 ", \
     2.69632693710332*x title " {/Times-Italic A}{/Symbol \242} = 2.69633 ", \
     2.5384699066375*x title " {/Times-Italic A}{/Symbol \242} = 2.53847 ", \
     1.27223718789601*x title " {/Times-Italic A}{/Symbol \242} = 1.27224 ", \
     2.76622165212169*x title " {/Times-Italic A}{/Symbol \242} = 2.76622 ", \
     2.2109717736501*x title " {/Times-Italic A}{/Symbol \242} = 2.21097 ", \
     1.53125978075883*x title " {/Times-Italic A}{/Symbol \242} = 1.53126 "
