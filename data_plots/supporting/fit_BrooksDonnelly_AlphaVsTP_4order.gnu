#filename: fit_BrooksDonnelly_AlphaVsTP_4order.gnu
reset
FIT_LIMIT = 1e-7
alph(x,y) = alph00      + alph10*x      + alph20*x**2      + alph30*x**3      + alph40*x**4      + alph01*y    + alph11*x*y    + alph21*x**2*y    + alph31*x**3*y    + alph41*x**4*y    + alph02*y**2 + alph12*x*y**2 + alph22*x**2*y**2 + alph32*x**3*y**2 + alph42*x**4*y**2 + alph03*y**3 + alph13*x*y**3 + alph23*x**2*y**3 + alph33*x**3*y**3 + alph43*x**4*y**3 + alph04*y**4 + alph14*x*y**4 + alph24*x**2*y**4 + alph34*x**3*y**4 + alph44*x**4*y**4
fit alph(x,y) 'BrooksDonnelly_AlphaVsTP.dat' using 1:($2*1.01325):4:(4) via 'fit_BrooksDonnelly_AlphaVsTP_4order_init.par'
# I first ran this fit using "via alph00,alph10,alph20,alph30,alph40, alph01,alph11,alph21,alph31,alph41, alph02,alph12,alph22,alph32,alph42, alph03,alph13,alph23,alph33,alph43, alph04,alph14,alph24,alph34,alph44" to get the first numbers for fit_BrooksDonnelly_AlphaVsTP_4order_init.par
update 'fit_BrooksDonnelly_AlphaVsTP_4order_init.par' 'fit_BrooksDonnelly_AlphaVsTP_4order_init.par'
set print 'fit_BrooksDonnelly_AlphaVsTP_4order.par'
print alph00,"\n", \
      alph10,"\n", \
      alph20,"\n", \
      alph30,"\n", \
      alph40,"\n", \
      alph01,"\n", \
      alph11,"\n", \
      alph21,"\n", \
      alph31,"\n", \
      alph41,"\n", \
      alph02,"\n", \
      alph12,"\n", \
      alph22,"\n", \
      alph32,"\n", \
      alph42,"\n", \
      alph03,"\n", \
      alph13,"\n", \
      alph23,"\n", \
      alph33,"\n", \
      alph43,"\n", \
      alph04,"\n", \
      alph14,"\n", \
      alph24,"\n", \
      alph34,"\n", \
      alph44
unset print
set print 'fit_BrooksDonnelly_AlphaVsTP_4order_init.par'
print "alph00 = ",alph00,"\n", \
      "alph10 = ",alph10,"\n", \
      "alph20 = ",alph20,"\n", \
      "alph30 = ",alph30,"\n", \
      "alph40 = ",alph40,"\n", \
      "alph01 = ",alph01,"\n", \
      "alph11 = ",alph11,"\n", \
      "alph21 = ",alph21,"\n", \
      "alph31 = ",alph31,"\n", \
      "alph41 = ",alph41,"\n", \
      "alph02 = ",alph02,"\n", \
      "alph12 = ",alph12,"\n", \
      "alph22 = ",alph22,"\n", \
      "alph32 = ",alph32,"\n", \
      "alph42 = ",alph42,"\n", \
      "alph03 = ",alph03,"\n", \
      "alph13 = ",alph13,"\n", \
      "alph23 = ",alph23,"\n", \
      "alph33 = ",alph33,"\n", \
      "alph43 = ",alph43,"\n", \
      "alph04 = ",alph04,"\n", \
      "alph14 = ",alph14,"\n", \
      "alph24 = ",alph24,"\n", \
      "alph34 = ",alph34,"\n", \
      "alph44 = ",alph44
unset print
