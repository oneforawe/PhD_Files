#filename: fit_Ahlers_TcVsP.gnu
fTc(x)  = T0 + T1*x + T2*x**2 + T3*x**3 + T4*x**4
fit fTc(x) 'fit_Ahlers_TcVsP.dat' using 2:1 via T0,T1,T2,T3,T4
update 'fit_Ahlers_TcVsP.par' 'fit_Ahlers_TcVsP.par'
set print 'fit_Ahlers_TcVsP.par'
print T0,"\n", \
      T1,"\n", \
      T2,"\n", \
      T3,"\n", \
      T4
unset print
