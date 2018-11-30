#filename: plot_3Dvlt_macro_Tv_Op12_lmax100_dl0.0001_RhosoRhoVsTv.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_3Dvlt_macro_Tv_Op12_lmax100_dl0.0001_RhosoRhoVsTv.ps'
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

plot '3Dvlt_macro_Tv_Op12_Cc0.20_lmax100_dl0.0001.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.20 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.25_lmax100_dl0.0001.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.25 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.30_lmax100_dl0.0001.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.30 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.35_lmax100_dl0.0001.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.35 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.40_lmax100_dl0.0001.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.40 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.45_lmax100_dl0.0001.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.45 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.50_lmax100_dl0.0001.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.50 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.55_lmax100_dl0.0001.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.55 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.60_lmax100_dl0.0001.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.60 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.65_lmax100_dl0.0001.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.65 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.70_lmax100_dl0.0001.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.70 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.75_lmax100_dl0.0001.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.75 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.80_lmax100_dl0.0001.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.80 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.85_lmax100_dl0.0001.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.85 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.90_lmax100_dl0.0001.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.90 ", \
     '3Dvlt_macro_Tv_Op12_Cc0.95_lmax100_dl0.0001.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.95 ", \
     '3Dvlt_macro_Tv_Op12_Cc1.00_lmax100_dl0.0001.dat' using 4:10 with points title " {/Times-Italic C}_c = 0.00 ", \
     '3Dvlt_macro_Tv_Op12_Cc1.05_lmax100_dl0.0001.dat' using 4:10 with points title " {/Times-Italic C}_c = 1.05 ", \
     '3Dvlt_macro_Tv_Op12_Cc1.10_lmax100_dl0.0001.dat' using 4:10 with points title " {/Times-Italic C}_c = 1.10 ", \
     '3Dvlt_macro_Tv_Op12_Cc1.15_lmax100_dl0.0001.dat' using 4:10 with points title " {/Times-Italic C}_c = 1.15 ", \
     '3Dvlt_macro_Tv_Op12_Cc1.20_lmax100_dl0.0001.dat' using 4:10 with points title " {/Times-Italic C}_c = 1.20 ", \
     0.947239578322491*x**0.67168851012355 title " 0.94724*{/Symbol-Oblique t}^{0.671689} ", \
     1.10530941978064*x**0.671688366345817 title " 1.10531*{/Symbol-Oblique t}^{0.671688} ", \
     1.24549829811706*x**0.671688210271479 title " 1.24550*{/Symbol-Oblique t}^{0.671688} ", \
     1.37150568643642*x**0.671688055569549 title " 1.37151*{/Symbol-Oblique t}^{0.671688} ", \
     1.48593675829359*x**0.671687901684964 title " 1.48594*{/Symbol-Oblique t}^{0.671688} ", \
     1.5907108814637*x**0.671687751408015 title " 1.59071*{/Symbol-Oblique t}^{0.671688} ", \
     1.68729117801053*x**0.671687606245146 title " 1.68729*{/Symbol-Oblique t}^{0.671688} ", \
     1.77682321475636*x**0.671687468175834 title " 1.77682*{/Symbol-Oblique t}^{0.671687} ", \
     1.86022334252808*x**0.671687334326688 title " 1.86022*{/Symbol-Oblique t}^{0.671687} ", \
     1.93823778746378*x**0.671687204315957 title " 1.93824*{/Symbol-Oblique t}^{0.671687} ", \
     2.01148336966919*x**0.67168708014834 title " 2.01148*{/Symbol-Oblique t}^{0.671687} ", \
     2.08047631748937*x**0.671686959954976 title " 2.08048*{/Symbol-Oblique t}^{0.671687} ", \
     2.14565353289856*x**0.671686844123217 title " 2.14565*{/Symbol-Oblique t}^{0.671687} ", \
     2.20738826512615*x**0.671686732295071 title " 2.20739*{/Symbol-Oblique t}^{0.671687} ", \
     2.26600213006192*x**0.671686625942894 title " 2.26600*{/Symbol-Oblique t}^{0.671687} ", \
     2.3217741117876*x**0.67168651989538 title " 2.32177*{/Symbol-Oblique t}^{0.671687} ", \
     2.37494832587502*x**0.671686422479029 title " 2.37495*{/Symbol-Oblique t}^{0.671686} ", \
     2.42573890619995*x**0.671686324584755 title " 2.42574*{/Symbol-Oblique t}^{0.671686} ", \
     2.47433543555489*x**0.671686230204483 title " 2.47434*{/Symbol-Oblique t}^{0.671686} ", \
     2.52090617002538*x**0.671686138309062 title " 2.52091*{/Symbol-Oblique t}^{0.671686} ", \
     2.56560133839719*x**0.671686051723933 title " 2.56560*{/Symbol-Oblique t}^{0.671686} "

# From the dl = 1e-05 version, we get essentially the same results:
#     0.947243000304917*x**0.671688772196827 title " 0.94724*{/Symbol-Oblique t}^{0.671689} ", \
#     1.10530945776666*x**0.671688387843576 title " 1.10531*{/Symbol-Oblique t}^{0.671688} ", \
#     1.24550557036177*x**0.671688696095295 title " 1.24551*{/Symbol-Oblique t}^{0.671689} ", \
#     1.37150898231895*x**0.67168823415248 title " 1.37151*{/Symbol-Oblique t}^{0.671688} ", \
#     1.48594109903934*x**0.671688106572027 title " 1.48594*{/Symbol-Oblique t}^{0.671688} ", \
#     1.59070922050506*x**0.671687641021204 title " 1.59071*{/Symbol-Oblique t}^{0.671688} ", \
#     1.68729874155136*x**0.671687976528042 title " 1.68730*{/Symbol-Oblique t}^{0.671688} ", \
#     1.77683113262539*x**0.67168785354224 title " 1.77683*{/Symbol-Oblique t}^{0.671688} ", \
#     1.86022734242609*x**0.671687481756715 title " 1.86023*{/Symbol-Oblique t}^{0.671687} ", \
#     1.93823267674893*x**0.671686980500563 title " 1.93823*{/Symbol-Oblique t}^{0.671687} ", \
#     2.0114880281104*x**0.671687254393285 title " 2.01149*{/Symbol-Oblique t}^{0.671687} ", \
#     2.08047585411995*x**0.671686950424154 title " 2.08048*{/Symbol-Oblique t}^{0.671687} ", \
#     2.14565589646719*x**0.671686932499539 title " 2.14566*{/Symbol-Oblique t}^{0.671687} ", \
#     2.20739411331556*x**0.671686947535253 title " 2.20739*{/Symbol-Oblique t}^{0.671687} ", \
#     2.26600232342521*x**0.671686641247804 title " 2.26600*{/Symbol-Oblique t}^{0.671687} ", \
#     2.32172939303938*x**0.671685029119356 title " 2.32173*{/Symbol-Oblique t}^{0.671685} ", \
#     2.37493297477569*x**0.67168592351716 title " 2.37493*{/Symbol-Oblique t}^{0.671686} ", \
#     2.42569847150951*x**0.671685004959187 title " 2.42570*{/Symbol-Oblique t}^{0.671685} ", \
#     2.47429605776997*x**0.671685009306275 title " 2.47430*{/Symbol-Oblique t}^{0.671685} ", \
#     2.52093440746836*x**0.671687028711283 title " 2.52093*{/Symbol-Oblique t}^{0.671687} ", \
#     2.56559745668519*x**0.671685913548401 title " 2.56560*{/Symbol-Oblique t}^{0.671686} "
