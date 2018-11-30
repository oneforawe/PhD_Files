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
     4.07565991200311*x**0.671695230092551 title " 4.07566*{/Symbol-Oblique t}^{0.671695} ", \
     3.90703402701179*x**0.671695562602823 title " 3.90703*{/Symbol-Oblique t}^{0.671696} ", \
     3.8347390528957*x**0.671695121237274 title " 3.83474*{/Symbol-Oblique t}^{0.671695} ", \
     3.81620959120735*x**0.671695139727103 title " 3.81621*{/Symbol-Oblique t}^{0.671695} ", \
     3.79750132788772*x**0.671695529383194 title " 3.79750*{/Symbol-Oblique t}^{0.671696} ", \
     3.77857875354631*x**0.671695471942009 title " 3.77858*{/Symbol-Oblique t}^{0.671695} ", \
     3.75948123848747*x**0.671696009798674 title " 3.75948*{/Symbol-Oblique t}^{0.671696} ", \
     3.7401256618963*x**0.671695313540152 title " 3.74013*{/Symbol-Oblique t}^{0.671695} ", \
     3.7206357142573*x**0.671696162549953 title " 3.72064*{/Symbol-Oblique t}^{0.671696} ", \
     3.70087720322309*x**0.671695697201906 title " 3.70088*{/Symbol-Oblique t}^{0.671696} ", \
     3.68090819087901*x**0.671695448768343 title " 3.68091*{/Symbol-Oblique t}^{0.671695} ", \
     3.66075251090048*x**0.671695711027683 title " 3.66075*{/Symbol-Oblique t}^{0.671696} ", \
     3.51310502308858*x**0.671694886434503 title " 3.51311*{/Symbol-Oblique t}^{0.671695} ", \
     3.28055067772942*x**0.671696243664281 title " 3.28055*{/Symbol-Oblique t}^{0.671696} ", \
     3.01739411021094*x**0.671691881516688 title " 3.01739*{/Symbol-Oblique t}^{0.671692} ", \
     2.68358012369547*x**0.668589211506349 title " 2.68358*{/Symbol-Oblique t}^{0.668589} ", \
     3.10567498000823*x**0.680609727414235 title " 3.10567*{/Symbol-Oblique t}^{0.680610} ", \
     2.68631283809922*x**0.668149060987026 title " 2.68631*{/Symbol-Oblique t}^{0.668149} ", \
     2.80948832702534*x**0.673471020057165 title " 2.80949*{/Symbol-Oblique t}^{0.673471} ", \
     2.74370324879451*x**0.676719159927542 title " 2.74370*{/Symbol-Oblique t}^{0.676719} ", \
     2.22982519596216*x**0.670980291098337 title " 2.22983*{/Symbol-Oblique t}^{0.670980} ", \
     1.80556253974676*x**0.681680353335755 title " 1.80556*{/Symbol-Oblique t}^{0.681680} "
