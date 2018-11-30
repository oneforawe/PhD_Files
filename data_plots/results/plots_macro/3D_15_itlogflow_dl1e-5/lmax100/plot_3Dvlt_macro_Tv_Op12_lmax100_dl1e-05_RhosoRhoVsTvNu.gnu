#filename: plot_3Dvlt_macro_Tv_Op12_lmax100_dl1e-05_RhosoRhoVsTvNu.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_3Dvlt_macro_Tv_Op12_lmax100_dl1e-05_RhosoRhoVsTvNu.ps'
set size 1.0,1.1

set title "^4He Superfluid Fraction {/Symbol-Oblique r}_s / {/Symbol-Oblique r} \n vs Exponentiated Temperature Variable {/Symbol-Oblique t}^{/Symbol-Oblique n}"

set ylabel "{/Symbol-Oblique r}_s / {/Symbol-Oblique r} (unitless)"
set xlabel "{/Symbol-Oblique t}^{/Symbol-Oblique n} = [{/Times-Italic 1 - T / T}_{/Symbol-Oblique l}{/Times-Italic (P)}]^{/Symbol-Oblique n} (unitless)"
#set xrange [0:3e-5]
#set yrange [0:1e-4]
set xrange [0:(3.6e-8)**0.671688352279845358]
set yrange [0:3e-5]
#set xrange [0:4.38e-4]
#set yrange [0:1e-3]
#set xrange [0:4e-5]
#set xrange [0:1e-5]
#set yrange [0:0.1]
set grid
set key outside right box width -16 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter 3Dvlt\\_macro.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.671688352279845358 "
set rmargin 35

plot '3Dvlt_macro_Tv_Op12_Cc0.20_lmax100_dl1e-05.dat' using ($4**0.671688352279845358):10 with points title " {/Times-Italic C}_c = 0.20 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.25_lmax100_dl1e-05.dat' using ($4**0.671688352279845358):10 with points title " {/Times-Italic C}_c = 0.25 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.30_lmax100_dl1e-05.dat' using ($4**0.671688352279845358):10 with points title " {/Times-Italic C}_c = 0.30 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.35_lmax100_dl1e-05.dat' using ($4**0.671688352279845358):10 with points title " {/Times-Italic C}_c = 0.35 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.40_lmax100_dl1e-05.dat' using ($4**0.671688352279845358):10 with points title " {/Times-Italic C}_c = 0.40 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.45_lmax100_dl1e-05.dat' using ($4**0.671688352279845358):10 with points title " {/Times-Italic C}_c = 0.45 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.50_lmax100_dl1e-05.dat' using ($4**0.671688352279845358):10 with points title " {/Times-Italic C}_c = 0.50 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.55_lmax100_dl1e-05.dat' using ($4**0.671688352279845358):10 with points title " {/Times-Italic C}_c = 0.55 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.60_lmax100_dl1e-05.dat' using ($4**0.671688352279845358):10 with points title " {/Times-Italic C}_c = 0.60 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.65_lmax100_dl1e-05.dat' using ($4**0.671688352279845358):10 with points title " {/Times-Italic C}_c = 0.65 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.70_lmax100_dl1e-05.dat' using ($4**0.671688352279845358):10 with points title " {/Times-Italic C}_c = 0.70 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.75_lmax100_dl1e-05.dat' using ($4**0.671688352279845358):10 with points title " {/Times-Italic C}_c = 0.75 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.80_lmax100_dl1e-05.dat' using ($4**0.671688352279845358):10 with points title " {/Times-Italic C}_c = 0.80 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.85_lmax100_dl1e-05.dat' using ($4**0.671688352279845358):10 with points title " {/Times-Italic C}_c = 0.85 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.90_lmax100_dl1e-05.dat' using ($4**0.671688352279845358):10 with points title " {/Times-Italic C}_c = 0.90 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.95_lmax100_dl1e-05.dat' using ($4**0.671688352279845358):10 with points title " {/Times-Italic C}_c = 0.95 ", \
     '3Dvlt_macro_Tv_Op12_Cc1.00_lmax100_dl1e-05.dat' using ($4**0.671688352279845358):10 with points title " {/Times-Italic C}_c = 1.00 ", \
     '3Dvlt_macro_Tv_Op12_Cc1.05_lmax100_dl1e-05.dat' using ($4**0.671688352279845358):10 with points title " {/Times-Italic C}_c = 1.05 ", \
     '3Dvlt_macro_Tv_Op12_Cc1.10_lmax100_dl1e-05.dat' using ($4**0.671688352279845358):10 with points title " {/Times-Italic C}_c = 1.10 ", \
     '3Dvlt_macro_Tv_Op12_Cc1.15_lmax100_dl1e-05.dat' using ($4**0.671688352279845358):10 with points title " {/Times-Italic C}_c = 1.15 ", \
     '3Dvlt_macro_Tv_Op12_Cc1.20_lmax100_dl1e-05.dat' using ($4**0.671688352279845358):10 with points title " {/Times-Italic C}_c = 1.20 ", \
     0.947236950765656*x title " {/Times-Italic A}{/Symbol \242} = 0.94724 ", \
     1.10530749672814*x  title " {/Times-Italic A}{/Symbol \242} = 1.10531 ", \
     1.24550003643335*x  title " {/Times-Italic A}{/Symbol \242} = 1.24550 ", \
     1.37151032891256*x  title " {/Times-Italic A}{/Symbol \242} = 1.37151 ", \
     1.48594573950362*x  title " {/Times-Italic A}{/Symbol \242} = 1.48595 ", \
     1.59072545003058*x  title " {/Times-Italic A}{/Symbol \242} = 1.59073 ", \
     1.68730948068671*x  title " {/Times-Italic A}{/Symbol \242} = 1.68731 ", \
     1.776844174785*x    title " {/Times-Italic A}{/Symbol \242} = 1.77684 ", \
     1.86024618307064*x  title " {/Times-Italic A}{/Symbol \242} = 1.86025 ", \
     1.93826779690228*x  title " {/Times-Italic A}{/Symbol \242} = 1.93827 ", \
     2.01151958934747*x  title " {/Times-Italic A}{/Symbol \242} = 2.01152 ", \
     2.08051399558705*x  title " {/Times-Italic A}{/Symbol \242} = 2.08051 ", \
     2.14569751672733*x  title " {/Times-Italic A}{/Symbol \242} = 2.14570 ", \
     2.20743515212332*x  title " {/Times-Italic A}{/Symbol \242} = 2.20744 ", \
     2.26605523199855*x  title " {/Times-Italic A}{/Symbol \242} = 2.26606 ", \
     2.32182850548036*x  title " {/Times-Italic A}{/Symbol \242} = 2.32183 ", \
     2.37501279045475*x  title " {/Times-Italic A}{/Symbol \242} = 2.37501 ", \
     2.42580659071287*x  title " {/Times-Italic A}{/Symbol \242} = 2.42581 ", \
     2.47440555007487*x  title " {/Times-Italic A}{/Symbol \242} = 2.47441 ", \
     2.5209797268657*x   title " {/Times-Italic A}{/Symbol \242} = 2.52098 ", \
     2.56568294465858*x  title " {/Times-Italic A}{/Symbol \242} = 2.56568 "
