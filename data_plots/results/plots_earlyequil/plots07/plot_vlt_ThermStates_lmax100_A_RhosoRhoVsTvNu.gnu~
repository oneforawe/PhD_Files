#filename: plot_vlt_ThermStates_lmax100_A_RhosoRhoVsTvNu.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_vlt_ThermStates_lmax100_A_RhosoRhoVsTvNu.ps'
set size 1.0,1.1

set title "^4He Superfluid Fraction {/Symbol-Oblique r}_s / {/Symbol-Oblique r} \n vs Exponentiated Temperature Variable {/Symbol-Oblique t}^{/Symbol-Oblique n}"

set ylabel "{/Symbol-Oblique r}_s / {/Symbol-Oblique r} (unitless)"
set xlabel "{/Symbol-Oblique t}^{/Symbol-Oblique n} = [{/Times-Italic 1 - T / T}_{/Symbol-Oblique l}{/Times-Italic (P)}]^{/Symbol-Oblique n} (unitless)"
set xrange [0:3e-5]
set yrange [0:1e-4]
#set xrange [0:(1e-5)**0.6716883522798452]
#set yrange [0:1e-3]
#set xrange [0:4.38e-4]
#set yrange [0:1e-3]
#set xrange [0:4e-5]
#set xrange [0:1e-5]
#set yrange [0:0.1]
set grid
set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
set rmargin 35

plot 'vlt_ThermStates_lmax100_A_Cc1.20.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 1.20 ", \
     'vlt_ThermStates_lmax100_A_Cc1.10.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 1.10 ", \
     'vlt_ThermStates_lmax100_A_Cc1.06.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 1.06 ", \
     'vlt_ThermStates_lmax100_A_Cc1.05.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 1.05 ", \
     'vlt_ThermStates_lmax100_A_Cc1.04.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 1.04 ", \
     'vlt_ThermStates_lmax100_A_Cc1.03.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 1.03 ", \
     'vlt_ThermStates_lmax100_A_Cc1.02.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 1.02 ", \
     'vlt_ThermStates_lmax100_A_Cc1.01.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 1.01 ", \
     'vlt_ThermStates_lmax100_A_Cc1.00.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 1.00 ", \
     'vlt_ThermStates_lmax100_A_Cc0.99.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.99 ", \
     'vlt_ThermStates_lmax100_A_Cc0.98.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.98 ", \
     'vlt_ThermStates_lmax100_A_Cc0.97.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.97 ", \
     'vlt_ThermStates_lmax100_A_Cc0.90.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.90 ", \
     'vlt_ThermStates_lmax100_A_Cc0.80.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.80 ", \
     'vlt_ThermStates_lmax100_A_Cc0.70.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.70 ", \
     'vlt_ThermStates_lmax100_A_Cc0.60.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.60 ", \
     'vlt_ThermStates_lmax100_A_Cc0.55.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.55 ", \
     'vlt_ThermStates_lmax100_A_Cc0.50.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.50 ", \
     'vlt_ThermStates_lmax100_A_Cc0.40.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.40 ", \
     'vlt_ThermStates_lmax100_A_Cc0.30.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.30 ", \
     'vlt_ThermStates_lmax100_A_Cc0.20.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.20 ", \
     'vlt_ThermStates_lmax100_A_Cc0.10.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.10 ", \
     2.55778101421577*x title " {/Times-Italic A}{/Symbol \242} = 2.55778 ", \
     2.46128272804751*x title " {/Times-Italic A}{/Symbol \242} = 2.46128 ", \
     2.42035278824969*x title " {/Times-Italic A}{/Symbol \242} = 2.42035 ", \
     2.40989605822985*x title " {/Times-Italic A}{/Symbol \242} = 2.40990 ", \
     2.39934102392742*x title " {/Times-Italic A}{/Symbol \242} = 2.39934 ", \
     2.38869914858582*x title " {/Times-Italic A}{/Symbol \242} = 2.38870 ", \
     2.37796254819872*x title " {/Times-Italic A}{/Symbol \242} = 2.37796 ", \
     2.36712413855473*x title " {/Times-Italic A}{/Symbol \242} = 2.36712 ", \
     2.35618912551054*x title " {/Times-Italic A}{/Symbol \242} = 2.35619 ", \
     2.34514983639514*x title " {/Times-Italic A}{/Symbol \242} = 2.34515 ", \
     2.3340112631923*x title " {/Times-Italic A}{/Symbol \242} = 2.33401 ", \
     2.32277194175202*x title " {/Times-Italic A}{/Symbol \242} = 2.32277 ", \
     2.24102355120801*x title " {/Times-Italic A}{/Symbol \242} = 2.24102 ", \
     2.11391733274349*x title " {/Times-Italic A}{/Symbol \242} = 2.11392 ", \
     1.97243780115601*x title " {/Times-Italic A}{/Symbol \242} = 1.97244 ", \
     1.81334876829543*x title " {/Times-Italic A}{/Symbol \242} = 1.81335 ", \
     1.72587287447208*x title " {/Times-Italic A}{/Symbol \242} = 1.72587 ", \
     1.63219525999867*x title " {/Times-Italic A}{/Symbol \242} = 1.63220 ", \
     1.42256867328662*x title " {/Times-Italic A}{/Symbol \242} = 1.42257 ", \
     1.17466952018548*x title " {/Times-Italic A}{/Symbol \242} = 1.17467 ", \
     0.872088507534735*x title " {/Times-Italic A}{/Symbol \242} = 0.87209 "
