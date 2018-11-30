#filename: plot_3Dvlt_macro_Tv_Op42_Kth4.147_lmax100_dl1e-05_RhosoRhoVsTv.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_3Dvlt_macro_Tv_Op42_Kth4.147_lmax100_dl1e-05_RhosoRhoVsTv.ps'
set size 1.3,1.1

set title "^4He Superfluid Fraction {/Symbol-Oblique r}_s / {/Symbol-Oblique r} \n vs Exponentiated Temperature Variable {/Symbol-Oblique t}^{/Symbol-Oblique n}"

set ylabel "{/Symbol-Oblique r}_s / {/Symbol-Oblique r} (unitless)"
set xlabel "{/Symbol-Oblique t} = [{/Times-Italic 1 - T / T}_{/Symbol-Oblique l}{/Times-Italic (P)}] (unitless)"
##set xrange [0:1e-5]
set xrange [0:5e-8]
#set xrange [0:4e-8]
#set yrange [0:3e-5]
#set yrange [0:0.9]
set grid
set key outside right box width -13 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter 3Dvlt\\_macro.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n}"
set rmargin 40

plot '3Dvlt_macro_Tv_Op42_Kth4.147_Cc0.20_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.20 ", \
     '3Dvlt_macro_Tv_Op42_Kth4.147_Cc0.25_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.25 ", \
     '3Dvlt_macro_Tv_Op42_Kth4.147_Cc0.30_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.30 ", \
     '3Dvlt_macro_Tv_Op42_Kth4.147_Cc0.35_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.35 ", \
     '3Dvlt_macro_Tv_Op42_Kth4.147_Cc0.40_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.40 ", \
     '3Dvlt_macro_Tv_Op42_Kth4.147_Cc0.45_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.45 ", \
     '3Dvlt_macro_Tv_Op42_Kth4.147_Cc0.50_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.50 ", \
     '3Dvlt_macro_Tv_Op42_Kth4.147_Cc0.55_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.55 ", \
     '3Dvlt_macro_Tv_Op42_Kth4.147_Cc0.60_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.60 ", \
     '3Dvlt_macro_Tv_Op42_Kth4.147_Cc0.65_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.65 ", \
     '3Dvlt_macro_Tv_Op42_Kth4.147_Cc0.70_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.70 ", \
     '3Dvlt_macro_Tv_Op42_Kth4.147_Cc0.75_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.75 ", \
     '3Dvlt_macro_Tv_Op42_Kth4.147_Cc0.80_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.80 ", \
     '3Dvlt_macro_Tv_Op42_Kth4.147_Cc0.85_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.85 ", \
     '3Dvlt_macro_Tv_Op42_Kth4.147_Cc0.90_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.90 ", \
     '3Dvlt_macro_Tv_Op42_Kth4.147_Cc0.95_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.95 ", \
     '3Dvlt_macro_Tv_Op42_Kth4.147_Cc1.00_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.00 ", \
     '3Dvlt_macro_Tv_Op42_Kth4.147_Cc1.05_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 1.05 ", \
     '3Dvlt_macro_Tv_Op42_Kth4.147_Cc1.10_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 1.10 ", \
     '3Dvlt_macro_Tv_Op42_Kth4.147_Cc1.15_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 1.15 ", \
     '3Dvlt_macro_Tv_Op42_Kth4.147_Cc1.20_lmax100_dl1e-05.dat' using 4:10 with points title " {/Times-Italic C}_c = 1.20 ", \
     0.940713259454108*x**0.67168851225439 title " 0.94071*{/Symbol-Oblique t}^{0.671689} ", \
     1.09769397733892*x**0.6716883646193 title " 1.09769*{/Symbol-Oblique t}^{0.671688} ", \
     1.2369169871178*x**0.671688209664756 title " 1.23692*{/Symbol-Oblique t}^{0.671688} ", \
     1.36205622327407*x**0.671688056188976 title " 1.36206*{/Symbol-Oblique t}^{0.671688} ", \
     1.47569885517289*x**0.671687900671721 title " 1.47570*{/Symbol-Oblique t}^{0.671688} ", \
     1.57975107617263*x**0.671687749271138 title " 1.57975*{/Symbol-Oblique t}^{0.671688} ", \
     1.67566594718066*x**0.6716876041925 title " 1.67567*{/Symbol-Oblique t}^{0.671688} ", \
     1.76458123350184*x**0.671687471550557 title " 1.76458*{/Symbol-Oblique t}^{0.671687} ", \
     1.84740667886366*x**0.671687334550354 title " 1.84741*{/Symbol-Oblique t}^{0.671687} ", \
     1.92488359275574*x**0.671687203512823 title " 1.92488*{/Symbol-Oblique t}^{0.671687} ", \
     1.99762448330129*x**0.671687077609323 title " 1.99762*{/Symbol-Oblique t}^{0.671687} ", \
     2.06614207963402*x**0.67168695743679 title " 2.06614*{/Symbol-Oblique t}^{0.671687} ", \
     2.13087031814994*x**0.671686845055122 title " 2.13087*{/Symbol-Oblique t}^{0.671687} ", \
     2.19217968406099*x**0.67168673223903 title " 2.19218*{/Symbol-Oblique t}^{0.671687} ", \
     2.2503895492381*x**0.67168662008479 title " 2.25039*{/Symbol-Oblique t}^{0.671687} ", \
     2.30577751844957*x**0.671686523255675 title " 2.30578*{/Symbol-Oblique t}^{0.671687} ", \
     2.3585851282426*x**0.671686417040016 title " 2.35859*{/Symbol-Oblique t}^{0.671686} ", \
     2.40902579870317*x**0.671686320380648 title " 2.40903*{/Symbol-Oblique t}^{0.671686} ", \
     2.45728771508391*x**0.671686233270368 title " 2.45729*{/Symbol-Oblique t}^{0.671686} ", \
     2.50353759679795*x**0.671686141666175 title " 2.50354*{/Symbol-Oblique t}^{0.671686} ", \
     2.54792474939171*x**0.671686052639668 title " 2.54792*{/Symbol-Oblique t}^{0.671686} "
