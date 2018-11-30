#filename: splot_BrooksDonnelly_AlphaVsTP_4order.gnu
reset
set terminal postscript color eps enhanced
set output 'splot_fit_BrooksDonnelly_AlphaVsTP_4order.ps'
#set size 1.3,1

set title "^4He Thermal Expansion Coefficient {/Symbol-Oblique a} \n vs Temperature {/Times-Italic T} and Pressure {/Times-Italic P} "

set zlabel "{/Symbol-Oblique a} (K^{-1})"
set xlabel "{/Times-Italic P} (bar)"
set ylabel "{/Times-Italic T} (K)"
#set logscale x
set xrange [0:2.5]
set yrange [0:30]
set zrange [-0.1:0.01]
#set grid
set key top box spacing 1.2 title " Fit of data from Brooks/Donnelly "
set lmargin 5

set pm3d
splot 'fit_BrooksDonnelly_AlphaVsTP_4order_3d.dat' with pm3d


#alph00 = -0.000438263632235238
#alph10 = 0.00335632456278826
#alph20 = -0.00703633174074553
#alph30 = 0.00735105423187598
#alph40 = -0.00309886444097797
#alph01 = -0.000152348446858034
#alph11 = 0.00125153930891024
#alph21 = -0.00262893683476524
#alph31 = 0.00172886076286486
#alph41 = -0.000446768599852708
#alph02 = 1.3538726063601e-05
#alph12 = -8.60154881727956e-05
#alph22 = 8.38683062146619e-05
#alph32 = 8.32314967470736e-05
#alph42 = -7.22072382346195e-05
#alph03 = -1.85603624862181e-06
#alph13 = 1.37049734997549e-05
#alph23 = -2.36237144589513e-05
#alph33 = 9.4945057952975e-06
#alph43 = 9.40860081678399e-07
#alph04 = 4.63778491416667e-08
#alph14 = -3.38327502527372e-07
#alph24 = 5.72279823076565e-07
#alph34 = -2.21840192950938e-07
#alph44 = -3.06162209919151e-08
#alph(x,y) = alph00      + alph10*x      + alph20*x**2      + alph30*x**3      + alph40*x**4      + alph01*y    + alph11*x*y    + alph21*x**2*y    + alph31*x**3*y    + alph41*x**4*y    + alph02*y**2 + alph12*x*y**2 + alph22*x**2*y**2 + alph32*x**3*y**2 + alph42*x**4*y**2 + alph03*y**3 + alph13*x*y**3 + alph23*x**2*y**3 + alph33*x**3*y**3 + alph43*x**4*y**3 + alph04*y**4 + alph14*x*y**4 + alph24*x**2*y**4 + alph34*x**3*y**4 + alph44*x**4*y**4

#set pm3d
#splot alph(x,y) with pm3d
