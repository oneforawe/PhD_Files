#filename: funkyfactorfit.gnu
#This is to get a function that will multiply Del/kB ("Delok") as a function of T/Tc to give the shape of curve that we'd like to see.
FIT_LIMIT = 1e-16
f(x) = a*exp(b*x-c)+0.99999
fit f(x) 'funkyfactorfit.dat' using 1:2 via a,b,c