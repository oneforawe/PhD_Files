#filename: plot_vlt_ThermStates_Dexp100_A_RhosoRhoVsTvAvgNu.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_vlt_ThermStates_Dexp100_A_RhosoRhoVsTvAvgNu.ps'
set size 1.3,1.1

set title "^4He Superfluid Fraction {/Symbol-Oblique r}_s / {/Symbol-Oblique r} \n vs Exponentiated Temperature Variable {/Symbol-Oblique t}^{/Symbol-Oblique n}"

set ylabel "{/Symbol-Oblique r}_s / {/Symbol-Oblique r} (unitless)"
set xlabel "{/Symbol-Oblique t}^{/Symbol-Oblique n} = [{/Times-Italic 1 - T / T}_{/Symbol-Oblique l}{/Times-Italic (P)}]^{0.67168386} (unitless)"
#set xrange [0:4e-8]
#set yrange [0:3e-5]
set grid
set key outside right box width -16 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n}"
set rmargin 21

plot 'vlt_ThermStates_Dexp100_A_Cc1.20.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 1.20 ", \
     'vlt_ThermStates_Dexp100_A_Cc1.10.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 1.10 ", \
     'vlt_ThermStates_Dexp100_A_Cc1.06.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 1.06 ", \
     'vlt_ThermStates_Dexp100_A_Cc1.05.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 1.05 ", \
     'vlt_ThermStates_Dexp100_A_Cc1.04.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 1.04 ", \
     'vlt_ThermStates_Dexp100_A_Cc1.03.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 1.03 ", \
     'vlt_ThermStates_Dexp100_A_Cc1.02.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 1.02 ", \
     'vlt_ThermStates_Dexp100_A_Cc1.01.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 1.01 ", \
     'vlt_ThermStates_Dexp100_A_Cc1.00.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 1.00 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.99.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 0.99 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.98.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 0.98 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.97.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 0.97 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.90.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 0.90 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.80.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 0.80 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.70.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 0.70 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.60.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 0.60 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.55.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 0.55 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.50.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 0.50 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.40.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 0.40 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.30.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 0.30 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.20.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 0.20 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.10.dat' using ($2**0.67168386):8 with points title " {/Times-Italic C}_c = 0.10 ", \
     2.54738231586686*x title " 2.54738*{/Symbol-Oblique t}^{0.67168386} ", \
     2.45675435257335*x title " 2.45675*{/Symbol-Oblique t}^{0.67168386} ", \
     2.41832542182998*x title " 2.41833*{/Symbol-Oblique t}^{0.67168386} ", \
     2.4085011279168*x title " 2.40850*{/Symbol-Oblique t}^{0.67168386} ", \
     2.39859476132646*x title " 2.39859*{/Symbol-Oblique t}^{0.67168386} ", \
     2.3886047859783*x title " 2.38860*{/Symbol-Oblique t}^{0.67168386} ", \
     2.37851916815073*x title " 2.37852*{/Symbol-Oblique t}^{0.67168386} ", \
     2.36834098874751*x title " 2.36834*{/Symbol-Oblique t}^{0.67168386} ", \
     2.35807140572823*x title " 2.35807*{/Symbol-Oblique t}^{0.67168386} ", \
     2.34770797688088*x title " 2.34771*{/Symbol-Oblique t}^{0.67168386} ", \
     2.3372489475025*x title " 2.33725*{/Symbol-Oblique t}^{0.67168386} ", \
     2.32668792237161*x title " 2.32669*{/Symbol-Oblique t}^{0.67168386} ", \
     2.24989315262881*x title " 2.24989*{/Symbol-Oblique t}^{0.67168386} ", \
     2.13039487265553*x title " 2.13039*{/Symbol-Oblique t}^{0.67168386} ", \
     1.99717365572264*x title " 1.99717*{/Symbol-Oblique t}^{0.67168386} ", \
     1.8469859215745*x title " 1.84699*{/Symbol-Oblique t}^{0.67168386} ", \
     1.76417756966128*x title " 1.76418*{/Symbol-Oblique t}^{0.67168386} ", \
     1.67528065881911*x title " 1.67528*{/Symbol-Oblique t}^{0.67168386} ", \
     1.47535542907093*x title " 1.47536*{/Symbol-Oblique t}^{0.67168386} ", \
     1.23662578348928*x title " 1.23663*{/Symbol-Oblique t}^{0.67168386} ", \
     0.940489589877759*x title " 0.94049*{/Symbol-Oblique t}^{0.67168386} ", \
     0.547473744433247*x title " 0.54747*{/Symbol-Oblique t}^{0.67168386} "
