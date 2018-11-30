#filename: fit_vlt_HeatCap_P_DK0_5e-07_i47_CpVsTv_PowForm.gnu
reset
#FIT_LIMIT = 1e-5
f(x) = c*x**0.0150649 + d
set print 'fit_vlt_HeatCap_P_DK0_5e-07_i47_CpVsTv_PowForm.par'
print "c","\t","d","\t","(in c_P = c*(1-T/T_lambda)**0.0150649 + d)","\n"
fit [1e-3:5e-6] f(x) 'vlt_HeatCap_P_00.050_DK0_5e-07_i47.dat' using 2:9 via c,d
print c,"\t",d
fit [1e-3:5e-6] f(x) 'vlt_HeatCap_P_01.646_DK0_5e-07_i47.dat' using 2:9 via c,d
print c,"\t",d
fit [1e-3:5e-6] f(x) 'vlt_HeatCap_P_07.328_DK0_5e-07_i47.dat' using 2:9 via c,d
print c,"\t",d
fit [1e-3:5e-6] f(x) 'vlt_HeatCap_P_15.031_DK0_5e-07_i47.dat' using 2:9 via c,d
print c,"\t",d
fit [1e-3:5e-6] f(x) 'vlt_HeatCap_P_18.180_DK0_5e-07_i47.dat' using 2:9 via c,d
print c,"\t",d
fit [1e-3:5e-6] f(x) 'vlt_HeatCap_P_22.533_DK0_5e-07_i47.dat' using 2:9 via c,d
print c,"\t",d
fit [1e-3:5e-6] f(x) 'vlt_HeatCap_P_25.868_DK0_5e-07_i47.dat' using 2:9 via c,d
print c,"\t",d
unset print
