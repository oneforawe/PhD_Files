#filename: HeDensityFits.gnu
#I guess this is meant for you to just copy and paste, adjusting the second number in "using n:m"
f(x) = a+b*x+c*x**2+d*x**3+e*x**4
fit f(x) 'HeliumDensitiesVsTempPress.dat' using 1:2 via a,b,c,d,e
