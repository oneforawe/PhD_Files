#filename: fit_Pdep_Cc.gnu
reset
fCc(x)  = C00 + C01*x + C02*x**2 + C03*x**3 + C04*x**4 + C05*x**5 + C06*x**6 + C07*x**7 + C08*x**8 + C09*x**9 + C10*x**10
fit fCc(x) 'Pdep_Cc_lmax100_dl0.001_Op2.dat' using 3:1 via C00,C01,C02,C03,C04,C05,C06,C07,C08,C09,C10
update 'fit_Pdep_Cc.par' 'fit_Pdep_Cc.par'
set print 'fit_Pdep_Cc.par'
print "#filename: fit_Pdep_Cc.par\n\n", \
      "Cc(P) fit parameters\n", \
      C00,"\n", \
      C01,"\n", \
      C02,"\n", \
      C03,"\n", \
      C04,"\n", \
      C05,"\n", \
      C06,"\n", \
      C07,"\n", \
      C08,"\n", \
      C09,"\n", \
      C10
unset print
set print 'fit_Pdep_Cc_results.dat'
print "filename: fit_Pdep_Cc_results.dat\n\n", \
      "Cc(P) fit parameters\n", \
      C00,"\n", \
      C01,"\n", \
      C02,"\n", \
      C03,"\n", \
      C04,"\n", \
      C05,"\n", \
      C06,"\n", \
      C07,"\n", \
      C08,"\n", \
      C09,"\n", \
      C10,"\n\n\n", \
      "If you want to plot the fit in gnuplot (with plot_fit_Pdep_Cc.gnu), use this:\n", \
      "#Defining the fit parameters and equation:\n", \
      "C00 = ",C00,"\n", \
      "C01 = ",C01,"\n", \
      "C02 = ",C02,"\n", \
      "C03 = ",C03,"\n", \
      "C04 = ",C04,"\n", \
      "C05 = ",C05,"\n", \
      "C06 = ",C06,"\n", \
      "C07 = ",C07,"\n", \
      "C08 = ",C08,"\n", \
      "C09 = ",C09,"\n", \
      "C10 = ",C10,"\n", \
      "fCc(x)  = C00 + C01*x + C02*x**2 + C03*x**3 + C04*x**4 + C05*x**5 + C06*x**6 + C07*x**7 + C08*x**8 + C09*x**9 + C10*x**10\n\n\n", \
      "So we have...\n", \
      "P (bar)\t","Cc(P)\n", \
      "0.050\t",fCc(0.050),"\n", \
      "1.646\t",fCc(1.646),"\n", \
      "7.328\t",fCc(7.328),"\n", \
      "15.031\t",fCc(15.031),"\n", \
      "18.180\t",fCc(18.180),"\n", \
      "22.533\t",fCc(22.533),"\n", \
      "25.868\t",fCc(25.868),"\n\n", \
      "P (bar=atm*1.01325)\tP (atm)\t","Cc(P)\n", \
      "0.0*1.01325\t0.0\t",fCc(0.0*1.01325),"\n", \
      "1.0*1.01325\t1.0\t",fCc(1.0*1.01325),"\n", \
      "5.0*1.01325\t5.0\t",fCc(5.0*1.01325),"\n", \
      "10.0*1.01325\t10.0\t",fCc(10.0*1.01325),"\n", \
      "15.0*1.01325\t15.0\t",fCc(15.0*1.01325),"\n", \
      "20.0*1.01325\t20.0\t",fCc(20.0*1.01325),"\n", \
      "24.0*1.01325\t24.0\t",fCc(24.0*1.01325),"\n\n", \
      "P (bar=atm*1.01325)\tP (atm)\t","Cc(P)\n", \
      "0.0*1.01325\t0.0\t",fCc(0.0*1.01325),"\n", \
      "2.5*1.01325\t2.5\t",fCc(2.5*1.01325),"\n", \
      "5.0*1.01325\t5.0\t",fCc(5.0*1.01325),"\n", \
      "7.5*1.01325\t7.5\t",fCc(7.5*1.01325),"\n", \
      "10.0*1.01325\t10.0\t",fCc(10.0*1.01325),"\n", \
      "12.5*1.01325\t12.5\t",fCc(12.5*1.01325),"\n", \
      "15.0*1.01325\t15.0\t",fCc(15.0*1.01325),"\n", \
      "17.5*1.01325\t17.5\t",fCc(17.5*1.01325),"\n", \
      "20.0*1.01325\t20.0\t",fCc(20.0*1.01325),"\n", \
      "22.5*1.01325\t22.5\t",fCc(22.5*1.01325),"\n", \
      "25.0*1.01325\t25.0\t",fCc(25.0*1.01325),"\n\n\n", \
      "Place the following in a vlt_CcK01Input.dat file and run vlt_K0cFind.c to find the K0c:\n\n", \
      "Cc\t","K01\n", \
      fCc(0.050),"\t",5,"\n", \
      fCc(1.646),"\t",5,"\n", \
      fCc(7.328),"\t",5,"\n", \
      fCc(15.031),"\t",5,"\n", \
      fCc(18.180),"\t",5,"\n", \
      fCc(22.533),"\t",5,"\n", \
      fCc(25.868),"\t",5,"\n\n", \
      "Cc\t","K01\n", \
      fCc(0.0*1.01325),"\t",5,"\n", \
      fCc(1.0*1.01325),"\t",5,"\n", \
      fCc(5.0*1.01325),"\t",5,"\n", \
      fCc(10.0*1.01325),"\t",5,"\n", \
      fCc(15.0*1.01325),"\t",5,"\n", \
      fCc(20.0*1.01325),"\t",5,"\n", \
      fCc(24.0*1.01325),"\t",5,"\n\n", \
      "Cc\t","K01\n", \
      fCc(0.0*1.01325),"\t",5,"\n", \
      fCc(2.5*1.01325),"\t",5,"\n", \
      fCc(5.0*1.01325),"\t",5,"\n", \
      fCc(7.5*1.01325),"\t",5,"\n", \
      fCc(10.0*1.01325),"\t",5,"\n", \
      fCc(12.5*1.01325),"\t",5,"\n", \
      fCc(15.0*1.01325),"\t",5,"\n", \
      fCc(17.5*1.01325),"\t",5,"\n", \
      fCc(20.0*1.01325),"\t",5,"\n", \
      fCc(22.5*1.01325),"\t",5,"\n", \
      fCc(25.0*1.01325),"\t",5,"\n"
unset print
