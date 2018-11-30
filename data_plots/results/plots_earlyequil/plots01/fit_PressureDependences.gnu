#filename: fit_PressureDependences.gnu
#FIT_LIMIT = 1e-16 #This fit doesn't seem to be sensitive to a change in FIT_LIMIT (maybe there's not enough data?)
fCc(x)  = C0 + C1*x + C2*x**2 + C3*x**3 + C4*x**4
fK0c(x) = K0 + K1*x + K2*x**2 + K3*x**3 + K4*x**4
# try adding more powers?
fit fCc(x)  'PressureDependences.dat' using ($3-2.46302)/-0.0281536:1 via C0,C1,C2,C3,C4
fit fK0c(x) 'PressureDependences.dat' using ($3-2.46302)/-0.0281536:2 via K0,K1,K2,K3,K4
update 'fit_PressureDependences.par' 'fit_PressureDependences.par'
set print 'fit_PressureDependences.par'
print "Parameters for Cc(P)\n", \
      C0,"\n", \
      C1,"\n", \
      C2,"\n", \
      C3,"\n", \
      C4,"\n\n", \
      "Parameters for K0c(P)\n", \
      K0,"\n", \
      K1,"\n", \
      K2,"\n", \
      K3,"\n", \
      K4
unset print
