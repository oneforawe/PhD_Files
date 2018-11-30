#filename: fit_3Dvlt_macro_cp_DK05e-07_i47_CpVsTv_PowForm2.gnu
reset
#FIT_LIMIT = 1e-5

alpha = 0.0150650568

f1(x) = c1*x**(-alpha) + d1
f2(x) = c2*x**(-alpha) + d2
f3(x) = c3*x**(-alpha) + d3
f4(x) = c4*x**(-alpha) + d4
f5(x) = c5*x**(-alpha) + d5
f6(x) = c6*x**(-alpha) + d6
f7(x) = c7*x**(-alpha) + d7

fit [1e-3:5e-6] f1(x) '3Dvlt_macro_cp_P00.050_DK05e-07_i47.dat' using 2:9 via c1,d1
fit [1e-3:5e-6] f2(x) '3Dvlt_macro_cp_P01.646_DK05e-07_i47.dat' using 2:9 via c2,d2
fit [1e-3:5e-6] f3(x) '3Dvlt_macro_cp_P07.328_DK05e-07_i47.dat' using 2:9 via c3,d3
fit [1e-3:5e-6] f4(x) '3Dvlt_macro_cp_P15.031_DK05e-07_i47.dat' using 2:9 via c4,d4
fit [1e-3:5e-6] f5(x) '3Dvlt_macro_cp_P18.180_DK05e-07_i47.dat' using 2:9 via c5,d5
fit [1e-3:5e-6] f6(x) '3Dvlt_macro_cp_P22.533_DK05e-07_i47.dat' using 2:9 via c6,d6
fit [1e-3:5e-6] f7(x) '3Dvlt_macro_cp_P25.868_DK05e-07_i47.dat' using 2:9 via c7,d7

set print 'fit_3Dvlt_macro_cp_DK05e-07_i47_CpVsTv_PowForm2.par'
print "filename: fit_3Dvlt_macro_cp_DK05e-07_i47_CpVsTv_PowForm2.par\n\n", \
      "c","\t","d","\t","(in c_P = c*(1-T/T_lambda)**(-alpha) + d)  where alpha = 0.0150650568","\n", \
      c1,"\t",d1,"\n", \
      c2,"\t",d2,"\n", \
      c3,"\t",d3,"\n", \
      c4,"\t",d4,"\n", \
      c5,"\t",d5,"\n", \
      c6,"\t",d6,"\n", \
      c7,"\t",d7,"\n\n\n\n", \
      "If you want to plot the fits in gnuplot (with plot_3Dvlt_macro_cp_PDK05e-07_i47_CpVsTv_AhlersCompare.gnu), use this:\n", \
      '     ',c1,'*x**(-alpha) + ',d1,' notitle, \',"\n", \
      '     ',c2,'*x**(-alpha) + ',d2,' notitle, \',"\n", \
      '     ',c3,'*x**(-alpha) + ',d3,' notitle, \',"\n", \
      '     ',c4,'*x**(-alpha) + ',d4,' notitle, \',"\n", \
      '     ',c5,'*x**(-alpha) + ',d5,' notitle, \',"\n", \
      '     ',c6,'*x**(-alpha) + ',d6,' notitle, \',"\n", \
      '     ',c7,'*x**(-alpha) + ',d7,' notitle'
unset print

a0 = 266.802325754247
a1 = 303.526950914442
a2 = 297.131603602319
a3 = 284.901964808795
a4 = 263.387041150176
a5 = 251.896132316713
a6 = 234.173581841199
a7 = 242.193065573482
b0 = -244.41096614849
b1 = -286.265877058042
b2 = -280.128048665732
b3 = -270.311080488692
b4 = -249.812444481508
b5 = -238.388959814274
b6 = -220.343521773895
b7 = -229.282735872321

set print 'HeatCapFits2.dat'
print "filename: HeatCapFits2.dat\n", \
      "You can copy and paste the following into HeatCapFits.ods:\n\n\n", \
      "SHIFTED POWER FORM FIT\n\n", \
      "\tc_P = A * (c_P)_vortex + (c_P)_ns\n", \
      "\t     (ns = nonsingular background part of specific heat)\n\n", \
      "\tc_P experimental data parameters               \t\t\t\t(c_P)_vortex calculation parameters\n", \
      "\t     ( c_P = a*(1-T/T_lambda)^(-alpha) + b )   \t\t\t\t     ( (c_P)_vortex = c*(1-T/T_lambda)^(-alpha) + d )\n\n\n", \
      "P\ta\tb\tc\td\t(c_P)_ns\tA\n", \
      "(bar)\t(J/(K.mol))\t(J/(K.mol))\t(J/(K.mol))\t(J/(K.mol))\t(J/(K.mol))\t(unitless)\n", \
      "0.05\t", a0,"\t",b0,"\t",c1,"\t",d1,"\t",(b0-a0*d1/c1),"\t",(a0/c1),"\n", \
      "1.65\t", a2,"\t",b2,"\t",c2,"\t",d2,"\t",(b2-a2*d2/c2),"\t",(a2/c2),"\n", \
      "7.33\t", a3,"\t",b3,"\t",c3,"\t",d3,"\t",(b3-a3*d3/c3),"\t",(a3/c3),"\n", \
      "15.03\t",a4,"\t",b4,"\t",c4,"\t",d4,"\t",(b4-a4*d4/c4),"\t",(a4/c4),"\n", \
      "18.18\t",a5,"\t",b5,"\t",c5,"\t",d5,"\t",(b5-a5*d5/c5),"\t",(a5/c5),"\n", \
      "22.53\t",a6,"\t",b6,"\t",c6,"\t",d6,"\t",(b6-a6*d6/c6),"\t",(a6/c6),"\n", \
      "25.87\t",a7,"\t",b7,"\t",c7,"\t",d7,"\t",(b7-a7*d7/c7),"\t",(a7/c7),"\n\n\n", \
      "Just the results:\n\n", \
      "(c_P)_ns\tA\n", \
      "(J/(K.mol))\t(unitless)\n", \
      (b0-a0*d1/c1),"\t",(a0/c1),"\n", \
      (b2-a2*d2/c2),"\t",(a2/c2),"\n", \
      (b3-a3*d3/c3),"\t",(a3/c3),"\n", \
      (b4-a4*d4/c4),"\t",(a4/c4),"\n", \
      (b5-a5*d5/c5),"\t",(a5/c5),"\n", \
      (b6-a6*d6/c6),"\t",(a6/c6),"\n", \
      (b7-a7*d7/c7),"\t",(a7/c7)
unset print

set print 'HeatCapFits2_4ThesisTable.dat'
print "filename: HeatCapFits2_4ThesisTable.dat\n", \
      "SHIFTED POWER FORM FIT\n\n", \
      "\tc_P = A * (c_P)_vortex + (c_P)_ns\n", \
      "\t     (ns = nonsingular background part of specific heat)\n\n", \
      "\tc_P experimental data parameters               \t\t\t\t(c_P)_vortex calculation parameters\n", \
      "\t     ( c_P = a*(1-T/T_lambda)^(-alpha) + b )   \t\t\t\t     ( (c_P)_vortex = c*(1-T/T_lambda)^(-alpha) + d )\n\n\n", \
      "P\ta\tb\tc\td\t(c_P)_ns\tA\n", \
      "(bar)\t(J/(K.mol))\t(J/(K.mol))\t(J/(K.mol))\t(J/(K.mol))\t(J/(K.mol))\t(unitless)\n", \
      " $0$&$05$ & $", sprintf("%4.1f",a0),"$ & $", sprintf("%4.1f",b0),"$ & $", sprintf("%5.1f",c1),"$ & $", sprintf("%5.1f",d1),"$ & $", sprintf("%4.2f",(b0-a0*d1/c1)),"$ & $", sprintf("%5.4f",(a0/c1)),"$ & $\n", \
      " $1$&$65$ & $", sprintf("%4.1f",a2),"$ & $", sprintf("%4.1f",b2),"$ & $", sprintf("%5.1f",c2),"$ & $", sprintf("%5.1f",d2),"$ & $", sprintf("%4.2f",(b2-a2*d2/c2)),"$ & $", sprintf("%5.4f",(a2/c2)),"$ & $\n", \
      " $7$&$33$ & $", sprintf("%4.1f",a3),"$ & $", sprintf("%4.1f",b3),"$ & $", sprintf("%5.1f",c3),"$ & $", sprintf("%5.1f",d3),"$ & $", sprintf("%4.2f",(b3-a3*d3/c3)),"$ & $", sprintf("%5.4f",(a3/c3)),"$ & $\n", \
      "$15$&$03$ & $", sprintf("%4.1f",a4),"$ & $", sprintf("%4.1f",b4),"$ & $", sprintf("%5.1f",c4),"$ & $", sprintf("%5.1f",d4),"$ & $", sprintf("%4.2f",(b4-a4*d4/c4)),"$ & $", sprintf("%5.4f",(a4/c4)),"$ & $\n", \
      "$18$&$18$ & $", sprintf("%4.1f",a5),"$ & $", sprintf("%4.1f",b5),"$ & $", sprintf("%5.1f",c5),"$ & $", sprintf("%5.1f",d5),"$ & $", sprintf("%4.2f",(b5-a5*d5/c5)),"$ & $", sprintf("%5.4f",(a5/c5)),"$ & $\n", \
      "$22$&$53$ & $", sprintf("%4.1f",a6),"$ & $", sprintf("%4.1f",b6),"$ & $", sprintf("%5.1f",c6),"$ & $", sprintf("%5.1f",d6),"$ & $", sprintf("%4.2f",(b6-a6*d6/c6)),"$ & $", sprintf("%5.4f",(a6/c6)),"$ & $\n", \
      "$25$&$87$ & $", sprintf("%4.1f",a7),"$ & $", sprintf("%4.1f",b7),"$ & $", sprintf("%5.1f",c7),"$ & $", sprintf("%5.1f",d7),"$ & $", sprintf("%4.2f",(b7-a7*d7/c7)),"$ & $", sprintf("%5.4f",(a7/c7)),"$ & $"
unset print
