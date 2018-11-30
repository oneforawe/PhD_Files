#filename: plot_vlt_ThermStates_Dexp100_A_RhosoRhoVsTvNu.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_vlt_ThermStates_Dexp100_A_RhosoRhoVsTvNu.ps'
set size 1.3,1.1

set title "^4He Superfluid Fraction {/Symbol-Oblique r}_s / {/Symbol-Oblique r} \n vs Exponentiated Temperature Variable {/Symbol-Oblique t}^{/Symbol-Oblique n}"

set ylabel "{/Symbol-Oblique r}_s / {/Symbol-Oblique r} (unitless)"
set xlabel "{/Symbol-Oblique t}^{/Symbol-Oblique n} = [{/Times-Italic 1 - T / T}_{/Symbol-Oblique l}({/Times-Italic P})]^{0.6716883} (unitless)"
#set xrange [0:1e-5]
#set yrange [0:0.01]
set grid
set key outside right box width -16 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n}"
set rmargin 24

plot 'vlt_ThermStates_Dexp100_A_beta_Cc1.20.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 1.20 ", \
     'vlt_ThermStates_Dexp100_A_beta_Cc1.10.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 1.10 ", \
     'vlt_ThermStates_Dexp100_A_beta_Cc1.06.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 1.06 ", \
     'vlt_ThermStates_Dexp100_A_beta_Cc1.05.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 1.05 ", \
     'vlt_ThermStates_Dexp100_A_beta_Cc1.04.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 1.04 ", \
     'vlt_ThermStates_Dexp100_A_beta_Cc1.03.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 1.03 ", \
     'vlt_ThermStates_Dexp100_A_beta_Cc1.02.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 1.02 ", \
     'vlt_ThermStates_Dexp100_A_beta_Cc1.01.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 1.01 ", \
     'vlt_ThermStates_Dexp100_A_beta_Cc1.00.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 1.00 ", \
     'vlt_ThermStates_Dexp100_A_beta_Cc0.99.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 0.99 ", \
     'vlt_ThermStates_Dexp100_A_beta_Cc0.98.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 0.98 ", \
     'vlt_ThermStates_Dexp100_A_beta_Cc0.97.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 0.97 ", \
     'vlt_ThermStates_Dexp100_A_beta_Cc0.90.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 0.90 ", \
     'vlt_ThermStates_Dexp100_A_beta_Cc0.80.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 0.80 ", \
     'vlt_ThermStates_Dexp100_A_beta_Cc0.70.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 0.70 ", \
     'vlt_ThermStates_Dexp100_A_beta_Cc0.60.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 0.60 ", \
     'vlt_ThermStates_Dexp100_A_beta_Cc0.55.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 0.55 ", \
     'vlt_ThermStates_Dexp100_A_beta_Cc0.50.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 0.50 ", \
     'vlt_ThermStates_Dexp100_A_beta_Cc0.40.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 0.40 ", \
     'vlt_ThermStates_Dexp100_A_beta_Cc0.30.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 0.30 ", \
     'vlt_ThermStates_Dexp100_A_beta_Cc0.20.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 0.20 ", \
     'vlt_ThermStates_Dexp100_A_beta_Cc0.10.dat' using ($2)**0.6716883:9 with points title " {/Times-Italic C}_c = 0.10 ", \
     2.54449351967991*x  title " {/Times-Italic A}{/Symbol \242} = 2.544 ", \
     2.45468784745952*x  title " {/Times-Italic A}{/Symbol \242} = 2.455 ", \
     2.41655181997362*x  title " {/Times-Italic A}{/Symbol \242} = 2.417 ", \
     2.40680717732104*x  title " {/Times-Italic A}{/Symbol \242} = 2.407 ", \
     2.39696907910705*x  title " {/Times-Italic A}{/Symbol \242} = 2.397 ", \
     2.38704286253215*x  title " {/Times-Italic A}{/Symbol \242} = 2.387 ", \
     2.37703142847134*x  title " {/Times-Italic A}{/Symbol \242} = 2.377 ", \
     2.36691782727221*x  title " {/Times-Italic A}{/Symbol \242} = 2.367 ", \
     2.35672123054231*x  title " {/Times-Italic A}{/Symbol \242} = 2.357 ", \
     2.34642163430972*x  title " {/Times-Italic A}{/Symbol \242} = 2.346 ", \
     2.33602284656519*x  title " {/Times-Italic A}{/Symbol \242} = 2.336 ", \
     2.32552885095095*x  title " {/Times-Italic A}{/Symbol \242} = 2.326 ", \
     2.24916096966919*x  title " {/Times-Italic A}{/Symbol \242} = 2.249 ", \
     2.13018349108229*x  title " {/Times-Italic A}{/Symbol \242} = 2.130 ", \
     1.99738344598232*x  title " {/Times-Italic A}{/Symbol \242} = 1.997 ", \
     1.8475207783541*x   title " {/Times-Italic A}{/Symbol \242} = 1.847 ", \
     1.76483851155854*x  title " {/Times-Italic A}{/Symbol \242} = 1.764 ", \
     1.67604460633597*x  title " {/Times-Italic A}{/Symbol \242} = 1.676 ", \
     1.47624816141727*x  title " {/Times-Italic A}{/Symbol \242} = 1.476 ", \
     1.23754017127143*x  title " {/Times-Italic A}{/Symbol \242} = 1.238 ", \
     0.941307318611049*x title " {/Times-Italic A}{/Symbol \242} = 0.941 ", \
     0.54804392773075*x  title " {/Times-Italic A}{/Symbol \242} = 0.548 "

