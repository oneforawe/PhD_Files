#filename: fit_BrooksDonnelly_HeRhoVsPT.gnu
reset
FIT_LIMIT = 1e-7
R(x,y) = R00      + R10*x      + R20*x**2      + R30*x**3      + R40*x**4      + R01*y    + R11*x*y    + R21*x**2*y    + R31*x**3*y    + R41*x**4*y    + R02*y**2 + R12*x*y**2 + R22*x**2*y**2 + R32*x**3*y**2 + R42*x**4*y**2 + R03*y**3 + R13*x*y**3 + R23*x**2*y**3 + R33*x**3*y**3 + R43*x**4*y**3 + R04*y**4 + R14*x*y**4 + R24*x**2*y**4 + R34*x**3*y**4 + R44*x**4*y**4
fit R(x,y) 'BrooksDonnelly_HeRhoVsTP_b.dat' using 1:($2*1.01325):4:(4) via 'fit_BrooksDonnelly_HeRhoVsPT_init.par'
# I first ran this fit using "via R00,R10,R20,R30,R40, R01,R11,R21,R31,R41, R02,R12,R22,R32,R42, R03,R13,R23,R33,R43, R04,R14,R24,R34,R44" to get the first numbers for fit_BrooksDonnelly_HeRhoVsPT_init.par
update 'fit_BrooksDonnelly_HeRhoVsPT_init.par' 'fit_BrooksDonnelly_HeRhoVsPT_init.par'
set print 'fit_BrooksDonnelly_HeRhoVsPT.par'
print R00,"\n", \
      R10,"\n", \
      R20,"\n", \
      R30,"\n", \
      R40,"\n", \
      R01,"\n", \
      R11,"\n", \
      R21,"\n", \
      R31,"\n", \
      R41,"\n", \
      R02,"\n", \
      R12,"\n", \
      R22,"\n", \
      R32,"\n", \
      R42,"\n", \
      R03,"\n", \
      R13,"\n", \
      R23,"\n", \
      R33,"\n", \
      R43,"\n", \
      R04,"\n", \
      R14,"\n", \
      R24,"\n", \
      R34,"\n", \
      R44
unset print
set print 'fit_BrooksDonnelly_HeRhoVsPT_init.par'
print "R00 = ",R00,"\n", \
      "R10 = ",R10,"\n", \
      "R20 = ",R20,"\n", \
      "R30 = ",R30,"\n", \
      "R40 = ",R40,"\n", \
      "R01 = ",R01,"\n", \
      "R11 = ",R11,"\n", \
      "R21 = ",R21,"\n", \
      "R31 = ",R31,"\n", \
      "R41 = ",R41,"\n", \
      "R02 = ",R02,"\n", \
      "R12 = ",R12,"\n", \
      "R22 = ",R22,"\n", \
      "R32 = ",R32,"\n", \
      "R42 = ",R42,"\n", \
      "R03 = ",R03,"\n", \
      "R13 = ",R13,"\n", \
      "R23 = ",R23,"\n", \
      "R33 = ",R33,"\n", \
      "R43 = ",R43,"\n", \
      "R04 = ",R04,"\n", \
      "R14 = ",R14,"\n", \
      "R24 = ",R24,"\n", \
      "R34 = ",R34,"\n", \
      "R44 = ",R44
unset print
