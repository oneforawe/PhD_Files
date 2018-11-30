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
     4.1195376212659*x title " {/Times-Italic A}{/Symbol \242} = 4.11954 ", \
     4.07732641390705*x title " {/Times-Italic A}{/Symbol \242} = 4.07733 ", \
     4.06004140102069*x title " {/Times-Italic A}{/Symbol \242} = 4.06004 ", \
     4.05568823360946*x title " {/Times-Italic A}{/Symbol \242} = 4.05569 ", \
     4.05131552540484*x title " {/Times-Italic A}{/Symbol \242} = 4.05132 ", \
     4.04692853127488*x title " {/Times-Italic A}{/Symbol \242} = 4.04693 ", \
     4.04252132469078*x title " {/Times-Italic A}{/Symbol \242} = 4.04252 ", \
     4.03810371611361*x title " {/Times-Italic A}{/Symbol \242} = 4.03810 ", \
     4.03367190413814*x title " {/Times-Italic A}{/Symbol \242} = 4.03367 ", \
     4.02922274343478*x title " {/Times-Italic A}{/Symbol \242} = 4.02922 ", \
     4.02475551741483*x title " {/Times-Italic A}{/Symbol \242} = 4.02476 ", \
     4.0202725305023*x title " {/Times-Italic A}{/Symbol \242} = 4.02027 ", \
     3.98845634328547*x title " {/Times-Italic A}{/Symbol \242} = 3.98846 ", \
     3.94156366416536*x title " {/Times-Italic A}{/Symbol \242} = 3.94156 ", \
     3.89285281945471*x title " {/Times-Italic A}{/Symbol \242} = 3.89285 ", \
     3.84216278924257*x title " {/Times-Italic A}{/Symbol \242} = 3.84216 ", \
     3.8160178570575*x title " {/Times-Italic A}{/Symbol \242} = 3.81602 ", \
     3.78931351575078*x title " {/Times-Italic A}{/Symbol \242} = 3.78931 ", \
     3.73409711084074*x title " {/Times-Italic A}{/Symbol \242} = 3.73410 ", \
     3.67627642204571*x title " {/Times-Italic A}{/Symbol \242} = 3.67628 ", \
     3.61557304304366*x title " {/Times-Italic A}{/Symbol \242} = 3.61557 ", \
     3.55164470040665*x title " {/Times-Italic A}{/Symbol \242} = 3.55164 "
