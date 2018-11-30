#filename: plot_3Dvlt_macro_Tv_lmax100_dl0.001_Op2_RhosoRhoVsTv.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_3Dvlt_macro_Tv_lmax100_dl0.001_Op2_RhosoRhoVsTv.ps'
set size 1.3,1.1

set title "^4He Superfluid Fraction {/Symbol-Oblique r}_s / {/Symbol-Oblique r} \n vs Exponentiated Temperature Variable {/Symbol-Oblique t}^{/Symbol-Oblique n}"

set ylabel "{/Symbol-Oblique r}_s / {/Symbol-Oblique r} (unitless)"
set xlabel "{/Symbol-Oblique t} = [{/Times-Italic 1 - T / T}_{/Symbol-Oblique l}{/Times-Italic (P)}] (unitless)"
#set xrange [0:1e-5]
set xrange [0:5e-8]
#set xrange [0:4e-8]
#set yrange [0:3e-5]
#set yrange [0:0.9]
set grid
set key outside right box width -13 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter 3Dvlt\\_macro.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n}"
set rmargin 40

plot '3Dvlt_macro_Tv_Cc1.20_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 1.20 ", \
     '3Dvlt_macro_Tv_Cc1.10_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 1.10 ", \
     '3Dvlt_macro_Tv_Cc1.06_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 1.06 ", \
     '3Dvlt_macro_Tv_Cc1.05_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 1.05 ", \
     '3Dvlt_macro_Tv_Cc1.04_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 1.04 ", \
     '3Dvlt_macro_Tv_Cc1.03_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 1.03 ", \
     '3Dvlt_macro_Tv_Cc1.02_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 1.02 ", \
     '3Dvlt_macro_Tv_Cc1.01_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 1.01 ", \
     '3Dvlt_macro_Tv_Cc1.00_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 1.00 ", \
     '3Dvlt_macro_Tv_Cc0.99_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.99 ", \
     '3Dvlt_macro_Tv_Cc0.98_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.98 ", \
     '3Dvlt_macro_Tv_Cc0.97_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.97 ", \
     '3Dvlt_macro_Tv_Cc0.90_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.90 ", \
     '3Dvlt_macro_Tv_Cc0.80_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.80 ", \
     '3Dvlt_macro_Tv_Cc0.70_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.70 ", \
     '3Dvlt_macro_Tv_Cc0.60_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.60 ", \
     '3Dvlt_macro_Tv_Cc0.55_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.55 ", \
     '3Dvlt_macro_Tv_Cc0.50_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.50 ", \
     '3Dvlt_macro_Tv_Cc0.40_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.40 ", \
     '3Dvlt_macro_Tv_Cc0.30_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.30 ", \
     '3Dvlt_macro_Tv_Cc0.20_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.20 ", \
     '3Dvlt_macro_Tv_Cc0.10_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.10 ", \
     4.01185985818162*x**0.671696129441382 title " 4.01186*{/Symbol-Oblique t}^{0.671696} ", \
     3.86670584158557*x**0.671695067146213 title " 3.86671*{/Symbol-Oblique t}^{0.671695} ", \
     3.80461954842125*x**0.671694890001726 title " 3.80462*{/Symbol-Oblique t}^{0.671695} ", \
     3.78873009392214*x**0.671695281236974 title " 3.78873*{/Symbol-Oblique t}^{0.671695} ", \
     3.77267590297189*x**0.671695718873332 title " 3.77268*{/Symbol-Oblique t}^{0.671696} ", \
     3.75644992356048*x**0.671696006455697 title " 3.75645*{/Symbol-Oblique t}^{0.671696} ", \
     3.73999894303144*x**0.671695180850566 title " 3.74000*{/Symbol-Oblique t}^{0.671695} ", \
     3.72345850214495*x**0.671696014455975 title " 3.72346*{/Symbol-Oblique t}^{0.671696} ", \
     3.70668595990894*x**0.671695658321957 title " 3.70669*{/Symbol-Oblique t}^{0.671696} ", \
     3.68977744428874*x**0.671696211180156 title " 3.68978*{/Symbol-Oblique t}^{0.671696} ", \
     3.67266888732972*x**0.671696319355596 title " 3.67267*{/Symbol-Oblique t}^{0.671696} ", \
     3.65530109798573*x**0.671695046145018 title " 3.65530*{/Symbol-Oblique t}^{0.671695} ", \
     3.52877781452682*x**0.671696978803974 title " 3.52878*{/Symbol-Oblique t}^{0.671697} ", \
     3.32898302352384*x**0.671696047415337 title " 3.32898*{/Symbol-Oblique t}^{0.671696} ", \
     3.10232886341127*x**0.671695384785829 title " 3.10233*{/Symbol-Oblique t}^{0.671695} ", \
     2.84192531013029*x**0.671695582555465 title " 2.84193*{/Symbol-Oblique t}^{0.671696} ", \
     2.69621108599739*x**0.671696383080435 title " 2.69621*{/Symbol-Oblique t}^{0.671696} ", \
     2.53800418023674*x**0.671693136422142 title " 2.53800*{/Symbol-Oblique t}^{0.671693} ", \
     2.42612793064319*x**0.679141480635883 title " 2.42613*{/Symbol-Oblique t}^{0.679141} ", \
     2.12573593255665*x**0.667347914275229 title " 2.12574*{/Symbol-Oblique t}^{0.667348} ", \
     2.03398155864283*x**0.670448816918122 title " 2.03398*{/Symbol-Oblique t}^{0.670449} ", \
     1.59429723678222*x**0.672276781122608 title " 1.59430*{/Symbol-Oblique t}^{0.672277} "
