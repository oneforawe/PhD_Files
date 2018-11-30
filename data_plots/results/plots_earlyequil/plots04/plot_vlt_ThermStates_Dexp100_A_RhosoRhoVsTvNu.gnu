#filename: plot_vlt_ThermStates_Dexp100_A_RhosoRhoVsTvNu.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_vlt_ThermStates_Dexp100_A_RhosoRhoVsTvNu.ps'
set size 1.0,1.1

set title "^4He Superfluid Fraction {/Symbol-Oblique r}_s / {/Symbol-Oblique r} \n vs Exponentiated Temperature Variable {/Symbol-Oblique t}^{/Symbol-Oblique n}"

set ylabel "{/Symbol-Oblique r}_s / {/Symbol-Oblique r} (unitless)"
set xlabel "{/Symbol-Oblique t}^{/Symbol-Oblique n} = [{/Times-Italic 1 - T / T}_{/Symbol-Oblique l}{/Times-Italic (P)}]^{/Symbol-Oblique n} (unitless)"
#set xrange [0:4e-5]
#set xrange [0:4.38e-4]
#set yrange [0:3e-5]
set grid
set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
set rmargin 35

plot 'vlt_ThermStates_Dexp100_A_Cc1.20.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 1.20 ", \
     'vlt_ThermStates_Dexp100_A_Cc1.10.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 1.10 ", \
     'vlt_ThermStates_Dexp100_A_Cc1.06.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 1.06 ", \
     'vlt_ThermStates_Dexp100_A_Cc1.05.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 1.05 ", \
     'vlt_ThermStates_Dexp100_A_Cc1.04.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 1.04 ", \
     'vlt_ThermStates_Dexp100_A_Cc1.03.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 1.03 ", \
     'vlt_ThermStates_Dexp100_A_Cc1.02.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 1.02 ", \
     'vlt_ThermStates_Dexp100_A_Cc1.01.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 1.01 ", \
     'vlt_ThermStates_Dexp100_A_Cc1.00.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 1.00 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.99.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.99 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.98.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.98 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.97.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.97 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.90.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.90 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.80.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.80 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.70.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.70 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.60.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.60 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.55.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.55 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.50.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.50 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.40.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.40 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.30.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.30 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.20.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.20 ", \
     'vlt_ThermStates_Dexp100_A_Cc0.10.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C}_c = 0.10 ", \
     2.5875700386783*x   title " {/Times-Italic A}{/Symbol \242} = 2.58757 ", \
     2.49965756709455*x  title " {/Times-Italic A}{/Symbol \242} = 2.49966 ", \
     2.4623598281966*x   title " {/Times-Italic A}{/Symbol \242} = 2.46236 ", \
     2.45282507965517*x  title " {/Times-Italic A}{/Symbol \242} = 2.45283 ", \
     2.44321107421416*x  title " {/Times-Italic A}{/Symbol \242} = 2.44321 ", \
     2.43350944874101*x  title " {/Times-Italic A}{/Symbol \242} = 2.43351 ", \
     2.42371822221908*x  title " {/Times-Italic A}{/Symbol \242} = 2.42372 ", \
     2.41383613494845*x  title " {/Times-Italic A}{/Symbol \242} = 2.41384 ", \
     2.40386957763097*x  title " {/Times-Italic A}{/Symbol \242} = 2.40387 ", \
     2.39380495498839*x  title " {/Times-Italic A}{/Symbol \242} = 2.39380 ", \
     2.38364622719623*x  title " {/Times-Italic A}{/Symbol \242} = 2.38365 ", \
     2.37338662493047*x  title " {/Times-Italic A}{/Symbol \242} = 2.37339 ", \
     2.29876088364832*x  title " {/Times-Italic A}{/Symbol \242} = 2.29876 ", \
     2.18247159792732*x  title " {/Times-Italic A}{/Symbol \242} = 2.18247 ", \
     2.05252360182846*x  title " {/Times-Italic A}{/Symbol \242} = 2.05252 ", \
     1.90551955468697*x  title " {/Times-Italic A}{/Symbol \242} = 1.90552 ", \
     1.82418680445609*x  title " {/Times-Italic A}{/Symbol \242} = 1.82419 ", \
     1.73661725597965*x  title " {/Times-Italic A}{/Symbol \242} = 1.73662 ", \
     1.53853793567444*x  title " {/Times-Italic A}{/Symbol \242} = 1.53854 ", \
     1.29951057127969*x  title " {/Times-Italic A}{/Symbol \242} = 1.29951 ", \
     0.998261085989499*x title " {/Times-Italic A}{/Symbol \242} = 0.99826 ", \
     0.588281851901256*x title " {/Times-Italic A}{/Symbol \242} = 0.58828 "
