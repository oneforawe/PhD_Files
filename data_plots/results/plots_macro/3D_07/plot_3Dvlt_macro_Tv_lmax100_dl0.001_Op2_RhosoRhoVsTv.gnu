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
     3.94679402460217*x**0.671695552132279 title " 3.94679*{/Symbol-Oblique t}^{0.671696} ", \
     3.84358155487515*x**0.671696262033354 title " 3.84358*{/Symbol-Oblique t}^{0.671696} ", \
     3.79984843604117*x**0.671694806382391 title " 3.79985*{/Symbol-Oblique t}^{0.671695} ", \
     3.7887405832111*x**0.671695523291551 title " 3.78874*{/Symbol-Oblique t}^{0.671696} ", \
     3.77750429784844*x**0.671695666510399 title " 3.77750*{/Symbol-Oblique t}^{0.671696} ", \
     3.76614765369102*x**0.671695150082356 title " 3.76615*{/Symbol-Oblique t}^{0.671695} ", \
     3.75469835671892*x**0.671694861707759 title " 3.75470*{/Symbol-Oblique t}^{0.671695} ", \
     3.74316944734918*x**0.671694913895825 title " 3.74317*{/Symbol-Oblique t}^{0.671695} ", \
     3.73156047139373*x**0.67169545693759 title " 3.73156*{/Symbol-Oblique t}^{0.671695} ", \
     3.71978039813032*x**0.671694516669903 title " 3.71978*{/Symbol-Oblique t}^{0.671695} ", \
     3.70797983450122*x**0.671695421175195 title " 3.70798*{/Symbol-Oblique t}^{0.671695} ", \
     3.69599517203528*x**0.671694619435854 title " 3.69600*{/Symbol-Oblique t}^{0.671695} ", \
     3.60917500392144*x**0.671694536923075 title " 3.60918*{/Symbol-Oblique t}^{0.671695} ", \
     3.47437184913104*x**0.671695805714103 title " 3.47437*{/Symbol-Oblique t}^{0.671696} ", \
     3.32407529163857*x**0.671695153020903 title " 3.32408*{/Symbol-Oblique t}^{0.671695} ", \
     3.15440790651999*x**0.671695963296786 title " 3.15441*{/Symbol-Oblique t}^{0.671696} ", \
     3.06054796266083*x**0.671696033859499 title " 3.06055*{/Symbol-Oblique t}^{0.671696} ", \
     2.9594084040773*x**0.671695280857085 title " 2.95941*{/Symbol-Oblique t}^{0.671695} ", \
     2.73019415643428*x**0.671695434037628 title " 2.73019*{/Symbol-Oblique t}^{0.671695} ", \
     2.45182827083413*x**0.671696427001284 title " 2.45183*{/Symbol-Oblique t}^{0.671696} ", \
     2.09609502324383*x**0.671695266986036 title " 2.09610*{/Symbol-Oblique t}^{0.671695} ", \
     1.59679610233322*x**0.671695411848993 title " 1.59680*{/Symbol-Oblique t}^{0.671695} "
