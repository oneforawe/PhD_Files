#filename: fit_3Dvlt_macro_cp_Tv_Op19_DK05e-07_lmax100_dl0.0001_CpVsTv_LogForm.gnu
reset
#FIT_LIMIT = 1e-5

# f(x) = c*log(x) + d
f1(x) = c1*log(x) + d1
f2(x) = c2*log(x) + d2
f3(x) = c3*log(x) + d3
f4(x) = c4*log(x) + d4
f5(x) = c5*log(x) + d5
f6(x) = c6*log(x) + d6
f7(x) = c7*log(x) + d7

fit [1e-3:5e-6] f1(x) '3Dvlt_macro_cp_Tv_Op19_P00.050_DK05e-07_lmax100_dl0.0001.dat' using 4:14 via c1,d1
fit [1e-3:5e-6] f2(x) '3Dvlt_macro_cp_Tv_Op19_P01.646_DK05e-07_lmax100_dl0.0001.dat' using 4:14 via c2,d2
fit [1e-3:5e-6] f3(x) '3Dvlt_macro_cp_Tv_Op19_P07.328_DK05e-07_lmax100_dl0.0001.dat' using 4:14 via c3,d3
fit [1e-3:5e-6] f4(x) '3Dvlt_macro_cp_Tv_Op19_P15.031_DK05e-07_lmax100_dl0.0001.dat' using 4:14 via c4,d4
fit [1e-3:5e-6] f5(x) '3Dvlt_macro_cp_Tv_Op19_P18.180_DK05e-07_lmax100_dl0.0001.dat' using 4:14 via c5,d5
fit [1e-3:5e-6] f6(x) '3Dvlt_macro_cp_Tv_Op19_P22.533_DK05e-07_lmax100_dl0.0001.dat' using 4:14 via c6,d6
fit [1e-3:5e-6] f7(x) '3Dvlt_macro_cp_Tv_Op19_P25.868_DK05e-07_lmax100_dl0.0001.dat' using 4:14 via c7,d7

set print 'fit_3Dvlt_macro_cp_Tv_Op19_DK05e-07_lmax100_dl0.0001_CpVsTv_LogForm.par'
print "filename: fit_3Dvlt_macro_cp_Tv_Op19_DK05e-07_lmax100_dl0.0001_CpVsTv_LogForm.par\n\n", \
      "c","\t","d","\t","(in c_P = c*log(1-T/T_lambda) + d)","\n", \
      c1,"\t",d1,"\n", \
      c2,"\t",d2,"\n", \
      c3,"\t",d3,"\n", \
      c4,"\t",d4,"\n", \
      c5,"\t",d5,"\n", \
      c6,"\t",d6,"\n", \
      c7,"\t",d7,"\n\n\n\n", \
      "If you want to plot the fits in gnuplot (with plot_3Dvlt_macro_cp_Tv_Op19_DK05e-07_lmax100_dl0.0001_CpVsTv_AhlersCompare.gnu), use this:\n", \
      '     ', c1, '*log(x) + ', d1, ' ls 0 title " Fits ", \',"\n", \
      '     ', c2, '*log(x) + ', d2, ' ls 0 notitle, \',"\n", \
      '     ', c3, '*log(x) + ', d3, ' ls 0 notitle, \',"\n", \
      '     ', c4, '*log(x) + ', d4, ' ls 0 notitle, \',"\n", \
      '     ', c5, '*log(x) + ', d5, ' ls 0 notitle, \',"\n", \
      '     ', c6, '*log(x) + ', d6, ' ls 0 notitle, \',"\n", \
      '     ', c7, '*log(x) + ', d7, ' ls 0 notitle'
unset print

# From fit_AhlersHeatCap_CpVsTv_LogForm3.par ...
a0 = -4.98533400515801
a1 = -5.25369492590319
a2 = -5.14640250160904
a3 = -4.80914613709055
a4 = -4.49297375716576
a5 = -4.29606180236529
a6 = -3.9914665285381
a7 = -4.12121735593117
b0 = 16.1955026214432
b1 = 14.2195710358391
b2 = 13.9223118097342
b3 = 12.63395581705
b4 = 11.4294158171749
b5 = 11.4479200157205
b6 = 11.9465372063032
b7 = 11.0395746010093

set print 'HeatCapFits_LogForm.dat'
print "filename: HeatCapFits_LogForm.dat\n", \
      "You can copy and paste the following into HeatCapFits.ods:\n\n\n", \
      "SHIFTED LOGARITHMIC FORM FIT\n\n", \
      "\tc_P = A * (c_P)_vortex + (c_P)_nl\n", \
      "\t     (nl = non-(vortex)-loop background part of specific heat)\n\n", \
      "\tc_P experimental data parameters        \t\t\t\t(c_P)_vortex calculation parameters\n", \
      "\t     ( c_P = a*ln(1-T/T_lambda) + b )   \t\t\t\t     ( (c_P)_vortex = c*ln(1-T/T_lambda) + d )\n\n\n", \
      "P\ta\tb\tc\td\t(c_P)_nl\tA\n", \
      "(bar)\t(J/(K.mol))\t(J/(K.mol))\t(J/(K.mol))\t(J/(K.mol))\t(J/(K.mol))\t(unitless)\n", \
      "0.05\t", a0,"\t",b0,"\t",c1,"\t",d1,"\t",(b0-a0*d1/c1),"\t",(a0/c1),"\n", \
      "1.65\t", a2,"\t",b2,"\t",c2,"\t",d2,"\t",(b2-a2*d2/c2),"\t",(a2/c2),"\n", \
      "7.33\t", a3,"\t",b3,"\t",c3,"\t",d3,"\t",(b3-a3*d3/c3),"\t",(a3/c3),"\n", \
      "15.03\t",a4,"\t",b4,"\t",c4,"\t",d4,"\t",(b4-a4*d4/c4),"\t",(a4/c4),"\n", \
      "18.18\t",a5,"\t",b5,"\t",c5,"\t",d5,"\t",(b5-a5*d5/c5),"\t",(a5/c5),"\n", \
      "22.53\t",a6,"\t",b6,"\t",c6,"\t",d6,"\t",(b6-a6*d6/c6),"\t",(a6/c6),"\n", \
      "25.87\t",a7,"\t",b7,"\t",c7,"\t",d7,"\t",(b7-a7*d7/c7),"\t",(a7/c7),"\n\n\n", \
      "Just the results:\n\n", \
      "(c_P)_nl\tA\n", \
      "(J/(K.mol))\t(unitless)\n", \
      (b0-a0*d1/c1),"\t",(a0/c1),"\n", \
      (b2-a2*d2/c2),"\t",(a2/c2),"\n", \
      (b3-a3*d3/c3),"\t",(a3/c3),"\n", \
      (b4-a4*d4/c4),"\t",(a4/c4),"\n", \
      (b5-a5*d5/c5),"\t",(a5/c5),"\n", \
      (b6-a6*d6/c6),"\t",(a6/c6),"\n", \
      (b7-a7*d7/c7),"\t",(a7/c7)
unset print

set print 'HeatCapFits_4ThesisTable.dat'
print "filename: HeatCapFits_4ThesisTable.dat\n", \
      "SHIFTED LOGARITHMIC FORM FIT\n\n", \
      "\tc_P = A * (c_P)_vortex + (c_P)_nl\n", \
      "\t     (nl = non-(vortex)-loop background part of specific heat)\n\n", \
      "\tc_P experimental data parameters        \t\t\t\t(c_P)_vortex calculation parameters\n", \
      "\t     ( c_P = a*ln(1-T/T_lambda) + b )   \t\t\t\t     ( (c_P)_vortex = c*ln(1-T/T_lambda) + d )\n\n\n", \
      "P\ta\tb\tc\td\t(c_P)_nl\tA\n", \
      "(bar)\t(J/(K.mol))\t(J/(K.mol))\t(J/(K.mol))\t(J/(K.mol))\t(J/(K.mol))\t(unitless)\n", \
      " $0$&$05$ & $", sprintf("%4.1f",a0),"$ & $", sprintf("%4.1f",b0),"$ & $", sprintf("%4.0f",c1),"$ & $", sprintf("%4.0f",d1),"$ & $", sprintf("%4.2f",(b0-a0*d1/c1)),"$ & $", sprintf("%5.4f",(a0/c1)),'$ & $$&$$  \\',"\n", \
      " $1$&$65$ & $", sprintf("%4.1f",a2),"$ & $", sprintf("%4.1f",b2),"$ & $", sprintf("%4.0f",c2),"$ & $", sprintf("%4.0f",d2),"$ & $", sprintf("%4.2f",(b2-a2*d2/c2)),"$ & $", sprintf("%5.4f",(a2/c2)),'$ & $$&$$  \\',"\n", \
      " $7$&$33$ & $", sprintf("%4.1f",a3),"$ & $", sprintf("%4.1f",b3),"$ & $", sprintf("%4.0f",c3),"$ & $", sprintf("%4.0f",d3),"$ & $", sprintf("%4.2f",(b3-a3*d3/c3)),"$ & $", sprintf("%5.4f",(a3/c3)),'$ & $$&$$  \\',"\n", \
      "$15$&$03$ & $", sprintf("%4.1f",a4),"$ & $", sprintf("%4.1f",b4),"$ & $", sprintf("%4.0f",c4),"$ & $", sprintf("%4.0f",d4),"$ & $", sprintf("%4.2f",(b4-a4*d4/c4)),"$ & $", sprintf("%5.4f",(a4/c4)),'$ & $$&$$  \\',"\n", \
      "$18$&$18$ & $", sprintf("%4.1f",a5),"$ & $", sprintf("%4.1f",b5),"$ & $", sprintf("%4.0f",c5),"$ & $", sprintf("%4.0f",d5),"$ & $", sprintf("%4.2f",(b5-a5*d5/c5)),"$ & $", sprintf("%5.4f",(a5/c5)),'$ & $$&$$  \\',"\n", \
      "$22$&$53$ & $", sprintf("%4.1f",a6),"$ & $", sprintf("%4.1f",b6),"$ & $", sprintf("%4.0f",c6),"$ & $", sprintf("%4.0f",d6),"$ & $", sprintf("%4.2f",(b6-a6*d6/c6)),"$ & $", sprintf("%5.4f",(a6/c6)),'$ & $$&$$  \\',"\n", \
      "$25$&$87$ & $", sprintf("%4.1f",a7),"$ & $", sprintf("%4.1f",b7),"$ & $", sprintf("%4.0f",c7),"$ & $", sprintf("%4.0f",d7),"$ & $", sprintf("%4.2f",(b7-a7*d7/c7)),"$ & $", sprintf("%5.4f",(a7/c7)),'$ & $$&$$  \bstrut\\'
unset print
