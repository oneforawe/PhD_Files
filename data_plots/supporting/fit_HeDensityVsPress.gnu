#filename: fit_HeDensityVsPress.gnu
f(x) = a + b*x + c*x**2 + d*x**3 + e*x**4
fit f(x) 'HeDensityTcVsPress.dat' using 1:2 via a,b,c,d,e
