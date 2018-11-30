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
     2.54748513421246*x**0.671687229602566 title " 2.54749*{/Symbol-Oblique t}^{0.671687} ", \
     2.45678567546662*x**0.671684942208145 title " 2.45679*{/Symbol-Oblique t}^{0.671685} ", \
     2.41843383975098*x**0.671687592834026 title " 2.41843*{/Symbol-Oblique t}^{0.671688} ", \
     2.40853171822087*x**0.671684918596759 title " 2.40853*{/Symbol-Oblique t}^{0.671685} ", \
     2.39862438659941*x**0.671684944026384 title " 2.39862*{/Symbol-Oblique t}^{0.671685} ", \
     2.38863273954357*x**0.671685094285988 title " 2.38863*{/Symbol-Oblique t}^{0.671685} ", \
     2.37862558683358*x**0.671687598201897 title " 2.37863*{/Symbol-Oblique t}^{0.671688} ", \
     2.36842644423992*x**0.671686867692743 title " 2.36843*{/Symbol-Oblique t}^{0.671687} ", \
     2.3581557161222*x**0.671686860006038 title " 2.35816*{/Symbol-Oblique t}^{0.671687} ", \
     2.34779956532926*x**0.671687124989894 title " 2.34780*{/Symbol-Oblique t}^{0.671687} ", \
     2.33736798757756*x**0.671688087118895 title " 2.33737*{/Symbol-Oblique t}^{0.671688} ", \
     2.32676055987661*x**0.671686453446591 title " 2.32676*{/Symbol-Oblique t}^{0.671686} ", \
     2.24996473305563*x**0.671686532756358 title " 2.24996*{/Symbol-Oblique t}^{0.671687} ", \
     2.1304777348989*x**0.671687091846755 title " 2.13048*{/Symbol-Oblique t}^{0.671687} ", \
     1.99723544537609*x**0.671686407351536 title " 1.99724*{/Symbol-Oblique t}^{0.671686} ", \
     1.84698240071057*x**0.671683747863081 title " 1.84698*{/Symbol-Oblique t}^{0.671684} ", \
     1.76412802862914*x**0.671681606626879 title " 1.76413*{/Symbol-Oblique t}^{0.671682} ", \
     1.67521344625374*x**0.671680606003616 title " 1.67521*{/Symbol-Oblique t}^{0.671681} ", \
     1.47525672944025*x**0.671678410467391 title " 1.47526*{/Symbol-Oblique t}^{0.671678} ", \
     1.23651403136313*x**0.671676494488984 title " 1.23651*{/Symbol-Oblique t}^{0.671676} ", \
     0.940449907525239*x**0.671680644172675 title " 0.94045*{/Symbol-Oblique t}^{0.671681} "

#     2.54748513421246*x**0.671687229602566 title " 2.54749*{/Symbol-Oblique t}^{0.671687} ", \
#     2.45678567546662*x**0.671684942208145 title " 2.45679*{/Symbol-Oblique t}^{0.671685} ", \
#     2.41843383975098*x**0.671687592834026 title " 2.41843*{/Symbol-Oblique t}^{0.671688} ", \
#     2.40853171822087*x**0.671684918596759 title " 2.40853*{/Symbol-Oblique t}^{0.671685} ", \
#     2.39862438659941*x**0.671684944026384 title " 2.39862*{/Symbol-Oblique t}^{0.671685} ", \
#     2.38863273954357*x**0.671685094285988 title " 2.38863*{/Symbol-Oblique t}^{0.671685} ", \
#     2.37862558683358*x**0.671687598201897 title " 2.37863*{/Symbol-Oblique t}^{0.671688} ", \
#     2.36842644423992*x**0.671686867692743 title " 2.36843*{/Symbol-Oblique t}^{0.671687} ", \
#     2.3581557161222*x**0.671686860006038 title " 2.35816*{/Symbol-Oblique t}^{0.671687} ", \
#     2.34779956532926*x**0.671687124989894 title " 2.34780*{/Symbol-Oblique t}^{0.671687} ", \
#     2.33736798757756*x**0.671688087118895 title " 2.33737*{/Symbol-Oblique t}^{0.671688} ", \
#     2.32676055987661*x**0.671686453446591 title " 2.32676*{/Symbol-Oblique t}^{0.671686} ", \
#     2.24996473305563*x**0.671686532756358 title " 2.24996*{/Symbol-Oblique t}^{0.671687} ", \
#     2.1304777348989*x**0.671687091846755 title " 2.13048*{/Symbol-Oblique t}^{0.671687} ", \
#     1.99723544537609*x**0.671686407351536 title " 1.99724*{/Symbol-Oblique t}^{0.671686} ", \
#     1.84698240071057*x**0.671683747863081 title " 1.84698*{/Symbol-Oblique t}^{0.671684} ", \
#     1.76412802862914*x**0.671681606626879 title " 1.76413*{/Symbol-Oblique t}^{0.671682} ", \
#     1.67521344625374*x**0.671680606003616 title " 1.67521*{/Symbol-Oblique t}^{0.671681} ", \
#     1.47525672944025*x**0.671678410467391 title " 1.47526*{/Symbol-Oblique t}^{0.671678} ", \
#     1.23651403136313*x**0.671676494488984 title " 1.23651*{/Symbol-Oblique t}^{0.671676} ", \
#     0.940449907525239*x**0.671680644172675 title " 0.94045*{/Symbol-Oblique t}^{0.671681} "

#     4.13047699281511*x**0.671694882328993 title " 4.13048*{/Symbol-Oblique t}^{0.671695} ", \
#     3.92296655069894*x**0.671694611192045 title " 3.92297*{/Symbol-Oblique t}^{0.671695} ", \
#     3.83436253314043*x**0.671695003394993 title " 3.83436*{/Symbol-Oblique t}^{0.671695} ", \
#     3.81155429029542*x**0.671692713478285 title " 3.81155*{/Symbol-Oblique t}^{0.671693} ", \
#     3.78872604741185*x**0.671694407787095 title " 3.78873*{/Symbol-Oblique t}^{0.671694} ", \
#     3.76556369106807*x**0.671693967364937 title " 3.76556*{/Symbol-Oblique t}^{0.671694} ", \
#     3.7421341776456*x**0.671692667873871 title " 3.74213*{/Symbol-Oblique t}^{0.671693} ", \
#     3.71858865409976*x**0.671693697133084 title " 3.71859*{/Symbol-Oblique t}^{0.671694} ", \
#     3.69477942033814*x**0.671694193030254 title " 3.69478*{/Symbol-Oblique t}^{0.671694} ", \
#     3.67065353822845*x**0.671692901614179 title " 3.67065*{/Symbol-Oblique t}^{0.671693} ", \
#     3.64637747230693*x**0.67169368558162 title " 3.64638*{/Symbol-Oblique t}^{0.671694} ", \
#     3.621816363425*x**0.671693473381176 title " 3.62182*{/Symbol-Oblique t}^{0.671693} ", \
#     3.49748164668918*x**0.670442178666503 title " 3.49748*{/Symbol-Oblique t}^{0.670442} ", \
#     3.60988193136314*x**0.672823966812658 title " 3.60988*{/Symbol-Oblique t}^{0.672824} ", \
#     3.52325963063991*x**0.670786452030357 title " 3.52326*{/Symbol-Oblique t}^{0.670786} ", \
#     3.69974712010191*x**0.676057425220422 title " 3.69975*{/Symbol-Oblique t}^{0.676057} ", \
#     3.40933962275802*x**0.670253539027993 title " 3.40934*{/Symbol-Oblique t}^{0.670254} ", \
#     3.38934966059891*x**0.671341437823637 title " 3.38935*{/Symbol-Oblique t}^{0.671341} ", \
#     3.16820819828331*x**0.670706982939665 title " 3.16821*{/Symbol-Oblique t}^{0.670707} ", \
#     2.87735380215567*x**0.671344714405332 title " 2.87735*{/Symbol-Oblique t}^{0.671345} ", \
#     2.53146641988492*x**0.676301437246539 title " 2.53147*{/Symbol-Oblique t}^{0.676301} "

#     '3Dvlt_macro_Tv_Cc0.10_lmax100_dl0.001_Op2.dat' using 2:6 with points title " {/Times-Italic C}_c = 0.10 ", \
