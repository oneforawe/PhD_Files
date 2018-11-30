#filename: DelokFitFake.gnu
reset
f(x,y) = a+b*y+c*exp(-(a+b*y)/x)*x/y+d*exp(-2*(a+b*y)/x)*y+(e+f*y+g*y**2)*exp(-3*(a+b*y)/x)
FIT_LIMIT = 1e-3
fit f(x,y) 'DelokFitFake.dat' using 1:2:3:(3) via 'parametersFake.par'
