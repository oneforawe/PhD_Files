#filename: plot_vlt_ThermStates_lmax100_A_RhosoRhoVsTvNu.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_vlt_ThermStates_lmax100_A_RhosoRhoVsTvNu.ps'
set size 1.0,1.1

set title "^4He Superfluid Fraction {/Symbol-Oblique r}_s / {/Symbol-Oblique r} \n vs Exponentiated Temperature Variable {/Symbol-Oblique t}^{/Symbol-Oblique n}"

set ylabel "{/Symbol-Oblique r}_s / {/Symbol-Oblique r} (unitless)"
set xlabel "{/Symbol-Oblique t}^{/Symbol-Oblique n} = [{/Times-Italic 1 - T / T}_{/Symbol-Oblique l}{/Times-Italic (P)}]^{/Symbol-Oblique n} (unitless)"
#set xrange [0:3e-5]
#set yrange [0:1e-4]
set xrange [0:(3.6e-8)**0.6716883522798452]
set yrange [0:3e-5]
#set xrange [0:4.38e-4]
#set yrange [0:1e-3]
#set xrange [0:4e-5]
#set xrange [0:1e-5]
#set yrange [0:0.1]
set grid
#set key outside right box width -20 height 0.5 spacing 1.2 title " Calculated from {/CM-Typewriter vlt\\_ThermStates.c} \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.6716883522798452"
set key outside right box width -3 height 0.5 spacing 1.2 title " VLT calculations \n Fit: {/Symbol-Oblique r}_s / {/Symbol-Oblique r} = {/Times-Italic A}{/Symbol \242}{/Symbol-Oblique t}^{/Symbol-Oblique n} \n {/Symbol-Oblique n} = 0.67168835 "
set rmargin 20

plot 'vlt_ThermStates_lmax100_A_Cc1.20.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C} = 1.20 ", \
     'vlt_ThermStates_lmax100_A_Cc1.10.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C} = 1.10 ", \
     'vlt_ThermStates_lmax100_A_Cc1.00.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C} = 1.00 ", \
     'vlt_ThermStates_lmax100_A_Cc0.90.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C} = 0.90 ", \
     'vlt_ThermStates_lmax100_A_Cc0.80.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C} = 0.80 ", \
     'vlt_ThermStates_lmax100_A_Cc0.70.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C} = 0.70 ", \
     'vlt_ThermStates_lmax100_A_Cc0.60.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C} = 0.60 ", \
     'vlt_ThermStates_lmax100_A_Cc0.50.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C} = 0.50 ", \
     'vlt_ThermStates_lmax100_A_Cc0.40.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C} = 0.40 ", \
     'vlt_ThermStates_lmax100_A_Cc0.30.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C} = 0.30 ", \
     'vlt_ThermStates_lmax100_A_Cc0.20.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C} = 0.20 ", \
     2.54752911319273*x ls 0 title " {/Times-Italic A}{/Symbol \242} = 2.54753 ", \
     2.45689950808939*x ls 0 title " {/Times-Italic A}{/Symbol \242} = 2.45690 ", \
     2.35821022192162*x ls 0 title " {/Times-Italic A}{/Symbol \242} = 2.35821 ", \
     2.25002122041447*x ls 0 title " {/Times-Italic A}{/Symbol \242} = 2.25002 ", \
     2.13051772135556*x ls 0 title " {/Times-Italic A}{/Symbol \242} = 2.13052 ", \
     1.99738510364482*x ls 0 title " {/Times-Italic A}{/Symbol \242} = 1.99739 ", \
     1.84752282696777*x ls 0 title " {/Times-Italic A}{/Symbol \242} = 1.84752 ", \
     1.67604671297598*x ls 0 title " {/Times-Italic A}{/Symbol \242} = 1.67605 ", \
     1.47624952349609*x ls 0 title " {/Times-Italic A}{/Symbol \242} = 1.47625 ", \
     1.23754156551565*x ls 0 title " {/Times-Italic A}{/Symbol \242} = 1.23754 ", \
     0.941308481751331*x ls 0 title " {/Times-Italic A}{/Symbol \242} = 0.94131 "

#     'vlt_ThermStates_lmax100_A_Cc1.03.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C} = 1.03 ", \
#     'vlt_ThermStates_lmax100_A_Cc1.02.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C} = 1.02 ", \
#     'vlt_ThermStates_lmax100_A_Cc1.01.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C} = 1.01 ", \
#     'vlt_ThermStates_lmax100_A_Cc0.99.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C} = 0.99 ", \
#     'vlt_ThermStates_lmax100_A_Cc0.98.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C} = 0.98 ", \
#     'vlt_ThermStates_lmax100_A_Cc0.97.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C} = 0.97 ", \
#     'vlt_ThermStates_lmax100_A_Cc0.55.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C} = 0.55 ", \
#     'vlt_ThermStates_lmax100_A_Cc0.10.dat' using ($2**0.6716883522798452):8 with points title " {/Times-Italic C} = 0.10 ", \
#     2.38873853395947*x title " {/Times-Italic A}{/Symbol \242} = 2.38874 ", \
#     2.37865611833976*x title " {/Times-Italic A}{/Symbol \242} = 2.37866 ", \
#     2.36847474728622*x title " {/Times-Italic A}{/Symbol \242} = 2.36847 ", \
#     2.34784163793297*x title " {/Times-Italic A}{/Symbol \242} = 2.34784 ", \
#     2.33738040150869*x title " {/Times-Italic A}{/Symbol \242} = 2.33738 ", \
#     2.32682244620833*x title " {/Times-Italic A}{/Symbol \242} = 2.32682 ", \
#     1.76484029862448*x title " {/Times-Italic A}{/Symbol \242} = 1.76484 ", \
