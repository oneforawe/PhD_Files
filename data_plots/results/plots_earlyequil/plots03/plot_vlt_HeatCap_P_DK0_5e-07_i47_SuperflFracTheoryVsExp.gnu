#filename: plot_vlt_HeatCap_P_DK0_5e-07_i47_SuperflFracTheoryVsExp.gnu
reset
set terminal postscript color eps enhanced
set output 'plot_vlt_HeatCap_P_DK0_5e-07_i47_SuperflFracTheoryVsExp.ps'

set title "^4He Superfluid Fraction {/Symbol-Oblique r}_s/{/Symbol-Oblique r} vs Unitless Temperature {/Times-Italic t} = {/Times-Italic T / T}_{/Symbol-Oblique l} \n at Various Pressures"

set ylabel "{/Symbol-Oblique r}_s/{/Symbol-Oblique r} (unitless)"
set xlabel "{/Times-Italic t} = {/Times-Italic T / T}_{/Symbol-Oblique l}{/Times-Italic (P)} (unitless)"
set xrange [0:1.1]
set yrange [0:1.1]
set grid
set key inside left bottom box width -8 spacing 1.2 title "Calculated from \n {/CM-Typewriter vlt\\_HeatCap\\_P.c}"
set rmargin 4

plot 'vlt_HeatCap_P_00.050_DK0_5e-07_i47.dat' using 1:8 with linespoints title " {/Times-Italic P} = {/Times-Italic P}_{sv} = 0.050 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.172 K) ", \
     'vlt_HeatCap_P_01.646_DK0_5e-07_i47.dat' using 1:8 with linespoints title " {/Times-Italic P} = 1.644 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.157 K) ", \
     'vlt_HeatCap_P_07.328_DK0_5e-07_i47.dat' using 1:8 with linespoints title " {/Times-Italic P} = 7.328 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.095 K) ", \
     'vlt_HeatCap_P_15.031_DK0_5e-07_i47.dat' using 1:8 with linespoints title " {/Times-Italic P} = 15.031 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.998 K) ", \
     'vlt_HeatCap_P_18.180_DK0_5e-07_i47.dat' using 1:8 with linespoints title " {/Times-Italic P} = 18.180 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.954 K) ", \
     'vlt_HeatCap_P_22.533_DK0_5e-07_i47.dat' using 1:8 with linespoints title " {/Times-Italic P} = 22.533 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.889 K) ", \
     'vlt_HeatCap_P_25.868_DK0_5e-07_i47.dat' using 1:8 with linespoints title " {/Times-Italic P} = 25.868 bar ({/Times-Italic T}_{/Symbol-Oblique l} = 2.836 K) "

#phonons: (2*pi**2*kB**4*T**4)/(45*hbar**3*c**5)
#rotons:  (2*sqrt[mu]*p0**4*exp[-Del/(kB*T)])/(3*(2*pi)**1.5*hbar**3*sqrt(kB*T))
# Pressure dependence of c, mu, p0, Del ?


J Meynard Phys Rev
Thermodynamic properties of He based on phonon/roton spectrum (in 70s?)