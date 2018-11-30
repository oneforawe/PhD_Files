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
     4.11864836047298*x**0.671692094695785 title " 4.11865*{/Symbol-Oblique t}^{0.671692} ", \
     4.07648251861469*x**0.671692558774066 title " 4.07648*{/Symbol-Oblique t}^{0.671693} ", \
     4.05916553275772*x**0.671691792733332 title " 4.05917*{/Symbol-Oblique t}^{0.671692} ", \
     4.0548177744434*x**0.671691923893306 title " 4.05482*{/Symbol-Oblique t}^{0.671692} ", \
     4.05047653513273*x**0.67169247082102 title " 4.05048*{/Symbol-Oblique t}^{0.671692} ", \
     4.04608026978766*x**0.671692295686718 title " 4.04608*{/Symbol-Oblique t}^{0.671692} ", \
     4.04168678511657*x**0.671692404720772 title " 4.04169*{/Symbol-Oblique t}^{0.671692} ", \
     4.03720140983813*x**0.671691116388689 title " 4.03720*{/Symbol-Oblique t}^{0.671691} ", \
     4.03278696459297*x**0.671691453185439 title " 4.03279*{/Symbol-Oblique t}^{0.671691} ", \
     4.02831127766855*x**0.671690865925681 title " 4.02831*{/Symbol-Oblique t}^{0.671691} ", \
     4.0239400922218*x**0.67169258873971 title " 4.02394*{/Symbol-Oblique t}^{0.671693} ", \
     4.01948560173514*x**0.671693024486273 title " 4.01949*{/Symbol-Oblique t}^{0.671693} ", \
     3.98765798233257*x**0.671692597703968 title " 3.98766*{/Symbol-Oblique t}^{0.671693} ", \
     3.94074325721897*x**0.671691869381833 title " 3.94074*{/Symbol-Oblique t}^{0.671692} ", \
     3.8920892839906*x**0.671692564825976 title " 3.89209*{/Symbol-Oblique t}^{0.671693} ", \
     3.84144915433117*x**0.671693089893118 title " 3.84145*{/Symbol-Oblique t}^{0.671693} ", \
     3.81526550557188*x**0.671692144292647 title " 3.81527*{/Symbol-Oblique t}^{0.671692} ", \
     3.78852233177533*x**0.671691146610901 title " 3.78852*{/Symbol-Oblique t}^{0.671691} ", \
     3.73342869860632*x**0.671693169814389 title " 3.73343*{/Symbol-Oblique t}^{0.671693} ", \
     3.67562295998318*x**0.671693088198355 title " 3.67562*{/Symbol-Oblique t}^{0.671693} ", \
     3.61489884172015*x**0.671692188454352 title " 3.61490*{/Symbol-Oblique t}^{0.671692} ", \
     3.55103038696389*x**0.671692805769015 title " 3.55103*{/Symbol-Oblique t}^{0.671693} "
