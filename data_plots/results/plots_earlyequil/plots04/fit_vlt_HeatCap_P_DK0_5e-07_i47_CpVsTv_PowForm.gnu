#filename: fit_vlt_HeatCap_P_DK0_5e-07_i47_CpVsTv_PowForm.gnu
reset
#FIT_LIMIT = 1e-5
f1(x) = c1*x**0.0150649 + d1
f2(x) = c2*x**0.0150649 + d2
f3(x) = c3*x**0.0150649 + d3
f4(x) = c4*x**0.0150649 + d4
f5(x) = c5*x**0.0150649 + d5
f6(x) = c6*x**0.0150649 + d6
f7(x) = c7*x**0.0150649 + d7
fit [1e-3:5e-6] f1(x) 'vlt_HeatCap_P_00.050_DK0_5e-07_i47.dat' using 2:9 via c1,d1
fit [1e-3:5e-6] f2(x) 'vlt_HeatCap_P_01.646_DK0_5e-07_i47.dat' using 2:9 via c2,d2
fit [1e-3:5e-6] f3(x) 'vlt_HeatCap_P_07.328_DK0_5e-07_i47.dat' using 2:9 via c3,d3
fit [1e-3:5e-6] f4(x) 'vlt_HeatCap_P_15.031_DK0_5e-07_i47.dat' using 2:9 via c4,d4
fit [1e-3:5e-6] f5(x) 'vlt_HeatCap_P_18.180_DK0_5e-07_i47.dat' using 2:9 via c5,d5
fit [1e-3:5e-6] f6(x) 'vlt_HeatCap_P_22.533_DK0_5e-07_i47.dat' using 2:9 via c6,d6
fit [1e-3:5e-6] f7(x) 'vlt_HeatCap_P_25.868_DK0_5e-07_i47.dat' using 2:9 via c7,d7
set print 'fit_vlt_HeatCap_P_DK0_5e-07_i47_CpVsTv_PowForm.par'
print "c","\t","d","\t","(in c_P = c*(1-T/T_lambda)**0.0150649 + d)","\n", \
      c1,"\t",d1,"\n", \
      c2,"\t",d2,"\n", \
      c3,"\t",d3,"\n", \
      c4,"\t",d4,"\n", \
      c5,"\t",d5,"\n", \
      c6,"\t",d6,"\n", \
      c7,"\t",d7,"\n\n\n\n", \
      '     ',c1,'*x**0.0150649 + ',d1,' notitle, \',"\n", \
      '     ',c2,'*x**0.0150649 + ',d2,' notitle, \',"\n", \
      '     ',c3,'*x**0.0150649 + ',d3,' notitle, \',"\n", \
      '     ',c4,'*x**0.0150649 + ',d4,' notitle, \',"\n", \
      '     ',c5,'*x**0.0150649 + ',d5,' notitle, \',"\n", \
      '     ',c6,'*x**0.0150649 + ',d6,' notitle, \',"\n", \
      '     ',c7,'*x**0.0150649 + ',d7,' notitle'
unset print
set print 'HeatCapFits.dat'
print "SHIFTED POWER FORM FIT\n\n", \
      "\tc_P = A * (c_P)_vortex + (c_P)_ns\n", \
      "\t     (ns = nonsingular background part of specific heat)\n\n", \
      "\tc_P experimental data parameters               \t\t\t\t(c_P)_vortex calculation parameters\n", \
      "\t     ( c_P = a*(1-T/T_lambda)^(-alpha) + b )   \t\t\t\t     ( (c_P)_vortex = c*(1-T/T_lambda)^(-alpha) + d )\n\n\n", \
      "P\ta\tb\tc\td\t(c_P)_ns\tA\n", \
      "(bar)\t(J/(K.mol))\t(J/(K.mol))\t(J/(K.mol))\t(J/(K.mol))\t(J/(K.mol))\t(unitless)\n", \
      "0.05\t-390.567788097653\t402.521480862073\t",c1,"\t",d1,"\t",(402.521480862073-(-390.567788097653)*d1/c1),"\t",(-390.567788097653/c1),"\n", \
      "1.65\t-382.152987948856\t393.974672534648\t",c2,"\t",d2,"\t",(393.974672534648-(-382.152987948856)*d1/c1),"\t",(-382.152987948856/c1),"\n", \
      "7.33\t-354.954707725270\t365.772775425809\t",c3,"\t",d3,"\t",(365.772775425809-(-354.954707725270)*d1/c1),"\t",(-354.954707725270/c1),"\n", \
      "15.03\t-337.535726960267\t346.648220580450\t",c4,"\t",d4,"\t",(346.648220580450-(-337.535726960267)*d1/c1),"\t",(-337.535726960267/c1),"\n", \
      "18.18\t-322.724551326792\t331.944462704508\t",c5,"\t",d5,"\t",(331.944462704508-(-322.724551326792)*d1/c1),"\t",(-322.724551326792/c1),"\n", \
      "22.53\t-299.582711936558\t309.500119848571\t",c6,"\t",d6,"\t",(309.500119848571-(-299.582711936558)*d1/c1),"\t",(-299.582711936558/c1),"\n", \
      "25.87\t-306.490059708448\t315.792579181738\t",c7,"\t",d7,"\t",(315.792579181738-(-306.490059708448)*d1/c1),"\t",(-306.490059708448/c1),"\n\n\n", \
      "Just the results:\n\n", \
      "(c_P)_ns\tA\n", \
      "(J/(K.mol))\t(unitless)\n", \
      (402.521480862073-(-390.567788097653)*d1/c1),"\t",(-390.567788097653/c1),"\n", \
      (393.974672534648-(-382.152987948856)*d1/c1),"\t",(-382.152987948856/c1),"\n", \
      (365.772775425809-(-354.954707725270)*d1/c1),"\t",(-354.954707725270/c1),"\n", \
      (346.648220580450-(-337.535726960267)*d1/c1),"\t",(-337.535726960267/c1),"\n", \
      (331.944462704508-(-322.724551326792)*d1/c1),"\t",(-322.724551326792/c1),"\n", \
      (309.500119848571-(-299.582711936558)*d1/c1),"\t",(-299.582711936558/c1),"\n", \
      (315.792579181738-(-306.490059708448)*d1/c1),"\t",(-306.490059708448/c1)
unset print
