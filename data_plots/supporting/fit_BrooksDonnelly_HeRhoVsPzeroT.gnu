#filename: fit_BrooksDonnelly_HeRhoVsPzeroT.gnu
R(x) = R0 + R1*x + R2*x**2 + R3*x**3 + R4*x**4
fit R(x) 'BrooksDonnelly_HeRhoVsTP.dat' using 1:2 via R0,R1,R2,R3,R4
update 'fit_BrooksDonnelly_HeRhoVsPzeroT.par' 'fit_BrooksDonnelly_HeRhoVsPzeroT.par'
set print 'fit_BrooksDonnelly_HeRhoVsPzeroT.par'
print R0,"\n", \
      R1,"\n", \
      R2,"\n", \
      R3,"\n", \
      R4
unset print
