#filename: plot_vlt_ThermStates_lmax100_A_RhosoRhoVsTvNu.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_vlt_ThermStates_lmax100_A_RhosoRhoVsTvNu.ps'
set size 1.0,1.1

set title "^4He Superfluid Fraction {/Symbol-Oblique r}_s / {/Symbol-Oblique r} \n vs Exponentiated Temperature Variable {/Symbol-Oblique t}^{/Symbol-Oblique n}"

set ylabel "{/Symbol-Oblique r}_s / {/Symbol-Oblique r} (unitless)"
set xlabel "{/Symbol-Oblique t}^{/Symbol-Oblique n} = [{/Times-Italic 1 - T / T}_{/Symbol-Oblique l}{/Times-Italic (P)}]^{/Symbol-Oblique n} (unitless)"
#set xrange [0:4.38e-4]
#set xrange [0:4e-5]
#set xrange [0:1e-5]
#set yrange [0:3e-5]
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
     12.534593710359*x title " {/Times-Italic A}{/Symbol \242} = 12.53459 ", \
     11.5913495183956*x title " {/Times-Italic A}{/Symbol \242} = 11.59135 ", \
     11.2040435600391*x title " {/Times-Italic A}{/Symbol \242} = 11.20404 ", \
     11.1063261581804*x title " {/Times-Italic A}{/Symbol \242} = 11.10633 ", \
     11.0082467678804*x title " {/Times-Italic A}{/Symbol \242} = 11.00825 ", \
     10.9097997876259*x title " {/Times-Italic A}{/Symbol \242} = 10.90980 ", \
     10.8110105909722*x title " {/Times-Italic A}{/Symbol \242} = 10.81101 ", \
     10.711853804364*x title " {/Times-Italic A}{/Symbol \242} = 10.71185 ", \
     10.6123426613863*x title " {/Times-Italic A}{/Symbol \242} = 10.61234 ", \
     10.512470466911*x title " {/Times-Italic A}{/Symbol \242} = 10.51247 ", \
     10.4122418855077*x title " {/Times-Italic A}{/Symbol \242} = 10.41224 ", \
     10.3116533462216*x title " {/Times-Italic A}{/Symbol \242} = 10.31165 ", \
     9.59766119990694*x title " {/Times-Italic A}{/Symbol \242} = 9.59766 ", \
     8.54829775380603*x title " {/Times-Italic A}{/Symbol \242} = 8.54830 ", \
     7.46657531467174*x title " {/Times-Italic A}{/Symbol \242} = 7.46658 ", \
     6.35684163136796*x title " {/Times-Italic A}{/Symbol \242} = 6.35684 ", \
     5.79359613610244*x title " {/Times-Italic A}{/Symbol \242} = 5.79360 ", \
     5.22638155805352*x title " {/Times-Italic A}{/Symbol \242} = 5.22638 ", \
     4.08687199318285*x title " {/Times-Italic A}{/Symbol \242} = 4.08687 ", \
     2.95641134473754*x title " {/Times-Italic A}{/Symbol \242} = 2.95641 ", \
     1.86224047542243*x title " {/Times-Italic A}{/Symbol \242} = 1.86224 ", \
     0.842254519808405*x title " {/Times-Italic A}{/Symbol \242} = 0.84225 "
