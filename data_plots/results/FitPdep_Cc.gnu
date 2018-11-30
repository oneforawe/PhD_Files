#filename: FitPdep_Cc.gnu
# Fit to find the pressure dependence of K0c
reset
f(x) = a+b*x+c*x**2+d*x**3+e*x**4
#FIT_LIMIT = 1e-3
fit f(x) 'RhosAmps.dat' using 4:1 via a,b,c,d,e
