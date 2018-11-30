#filename: fit_3Dvlt_macro_Tv_lmax100_dl0.001_Op2_RhosoRhoVsTvNu.gnu
f01(x) = A01*x**0.671696337558
f02(x) = A02*x**0.671696337558
f03(x) = A03*x**0.671696337558
f04(x) = A04*x**0.671696337558
f05(x) = A05*x**0.671696337558
f06(x) = A06*x**0.671696337558
f07(x) = A07*x**0.671696337558
f08(x) = A08*x**0.671696337558
f09(x) = A09*x**0.671696337558
f10(x) = A10*x**0.671696337558
f11(x) = A11*x**0.671696337558
f12(x) = A12*x**0.671696337558
f13(x) = A13*x**0.671696337558
f14(x) = A14*x**0.671696337558
f15(x) = A15*x**0.671696337558
f16(x) = A16*x**0.671696337558
f17(x) = A17*x**0.671696337558
f18(x) = A18*x**0.671696337558
f19(x) = A19*x**0.671696337558
f20(x) = A20*x**0.671696337558
f21(x) = A21*x**0.671696337558
fit [0:3.6e-8] f01(x) '3Dvlt_macro_Tv_Cc1.20_lmax100_dl0.001_Op2.dat' using 2:6 via A01
fit [0:3.6e-8] f02(x) '3Dvlt_macro_Tv_Cc1.10_lmax100_dl0.001_Op2.dat' using 2:6 via A02
fit [0:3.6e-8] f03(x) '3Dvlt_macro_Tv_Cc1.06_lmax100_dl0.001_Op2.dat' using 2:6 via A03
fit [0:3.6e-8] f04(x) '3Dvlt_macro_Tv_Cc1.05_lmax100_dl0.001_Op2.dat' using 2:6 via A04
fit [0:3.6e-8] f05(x) '3Dvlt_macro_Tv_Cc1.04_lmax100_dl0.001_Op2.dat' using 2:6 via A05
fit [0:3.6e-8] f06(x) '3Dvlt_macro_Tv_Cc1.03_lmax100_dl0.001_Op2.dat' using 2:6 via A06
fit [0:3.6e-8] f07(x) '3Dvlt_macro_Tv_Cc1.02_lmax100_dl0.001_Op2.dat' using 2:6 via A07
fit [0:3.6e-8] f08(x) '3Dvlt_macro_Tv_Cc1.01_lmax100_dl0.001_Op2.dat' using 2:6 via A08
fit [0:3.6e-8] f09(x) '3Dvlt_macro_Tv_Cc1.00_lmax100_dl0.001_Op2.dat' using 2:6 via A09
fit [0:3.6e-8] f10(x) '3Dvlt_macro_Tv_Cc0.99_lmax100_dl0.001_Op2.dat' using 2:6 via A10
fit [0:3.6e-8] f11(x) '3Dvlt_macro_Tv_Cc0.98_lmax100_dl0.001_Op2.dat' using 2:6 via A11
fit [0:3.6e-8] f12(x) '3Dvlt_macro_Tv_Cc0.97_lmax100_dl0.001_Op2.dat' using 2:6 via A12
fit [0:3.6e-8] f13(x) '3Dvlt_macro_Tv_Cc0.90_lmax100_dl0.001_Op2.dat' using 2:6 via A13
fit [0:3.6e-8] f14(x) '3Dvlt_macro_Tv_Cc0.80_lmax100_dl0.001_Op2.dat' using 2:6 via A14
fit [0:3.6e-8] f15(x) '3Dvlt_macro_Tv_Cc0.70_lmax100_dl0.001_Op2.dat' using 2:6 via A15
fit [0:3.6e-8] f16(x) '3Dvlt_macro_Tv_Cc0.60_lmax100_dl0.001_Op2.dat' using 2:6 via A16
fit [0:3.6e-8] f17(x) '3Dvlt_macro_Tv_Cc0.55_lmax100_dl0.001_Op2.dat' using 2:6 via A17
fit [0:3.6e-8] f18(x) '3Dvlt_macro_Tv_Cc0.50_lmax100_dl0.001_Op2.dat' using 2:6 via A18
fit [0:3.6e-8] f19(x) '3Dvlt_macro_Tv_Cc0.40_lmax100_dl0.001_Op2.dat' using 2:6 via A19
fit [0:3.6e-8] f20(x) '3Dvlt_macro_Tv_Cc0.30_lmax100_dl0.001_Op2.dat' using 2:6 via A20
fit [0:3.6e-8] f21(x) '3Dvlt_macro_Tv_Cc0.20_lmax100_dl0.001_Op2.dat' using 2:6 via A21
update 'fit_3Dvlt_macro_Tv_lmax100_dl0.001_Op2_RhosoRhoVsTvNu.par' 'fit_3Dvlt_macro_Tv_lmax100_dl0.001_Op2_RhosoRhoVsTvNu.par'
set print 'fit_3Dvlt_macro_Tv_lmax100_dl0.001_Op2_RhosoRhoVsTvNu.par'
print "A' ( for  Rhos/Rho = A'*(1-T/Tc)^nu  at different Cc's / K0c's / Pressures )\n", \
      A01,"\n", \
      A02,"\n", \
      A03,"\n", \
      A04,"\n", \
      A05,"\n", \
      A06,"\n", \
      A07,"\n", \
      A08,"\n", \
      A09,"\n", \
      A10,"\n", \
      A11,"\n", \
      A12,"\n", \
      A13,"\n", \
      A14,"\n", \
      A15,"\n", \
      A16,"\n", \
      A17,"\n", \
      A18,"\n", \
      A19,"\n", \
      A20,"\n", \
      A21,"\n\n\n", \
      "Put the following in plot_3Dvlt_macro_Tv_lmax100_dl0.001_Op2_RhosoRhoVsTvNu.gnu\n\n", \
      '     ',A01,'*x title " {/Times-Italic A}{/Symbol \242} = ',sprintf("%0.5f",A01),' ", \',"\n", \
      '     ',A02,'*x title " {/Times-Italic A}{/Symbol \242} = ',sprintf("%0.5f",A02),' ", \',"\n", \
      '     ',A03,'*x title " {/Times-Italic A}{/Symbol \242} = ',sprintf("%0.5f",A03),' ", \',"\n", \
      '     ',A04,'*x title " {/Times-Italic A}{/Symbol \242} = ',sprintf("%0.5f",A04),' ", \',"\n", \
      '     ',A05,'*x title " {/Times-Italic A}{/Symbol \242} = ',sprintf("%0.5f",A05),' ", \',"\n", \
      '     ',A06,'*x title " {/Times-Italic A}{/Symbol \242} = ',sprintf("%0.5f",A06),' ", \',"\n", \
      '     ',A07,'*x title " {/Times-Italic A}{/Symbol \242} = ',sprintf("%0.5f",A07),' ", \',"\n", \
      '     ',A08,'*x title " {/Times-Italic A}{/Symbol \242} = ',sprintf("%0.5f",A08),' ", \',"\n", \
      '     ',A09,'*x title " {/Times-Italic A}{/Symbol \242} = ',sprintf("%0.5f",A09),' ", \',"\n", \
      '     ',A10,'*x title " {/Times-Italic A}{/Symbol \242} = ',sprintf("%0.5f",A10),' ", \',"\n", \
      '     ',A11,'*x title " {/Times-Italic A}{/Symbol \242} = ',sprintf("%0.5f",A11),' ", \',"\n", \
      '     ',A12,'*x title " {/Times-Italic A}{/Symbol \242} = ',sprintf("%0.5f",A12),' ", \',"\n", \
      '     ',A13,'*x title " {/Times-Italic A}{/Symbol \242} = ',sprintf("%0.5f",A13),' ", \',"\n", \
      '     ',A14,'*x title " {/Times-Italic A}{/Symbol \242} = ',sprintf("%0.5f",A14),' ", \',"\n", \
      '     ',A15,'*x title " {/Times-Italic A}{/Symbol \242} = ',sprintf("%0.5f",A15),' ", \',"\n", \
      '     ',A16,'*x title " {/Times-Italic A}{/Symbol \242} = ',sprintf("%0.5f",A16),' ", \',"\n", \
      '     ',A17,'*x title " {/Times-Italic A}{/Symbol \242} = ',sprintf("%0.5f",A17),' ", \',"\n", \
      '     ',A18,'*x title " {/Times-Italic A}{/Symbol \242} = ',sprintf("%0.5f",A18),' ", \',"\n", \
      '     ',A19,'*x title " {/Times-Italic A}{/Symbol \242} = ',sprintf("%0.5f",A19),' ", \',"\n", \
      '     ',A20,'*x title " {/Times-Italic A}{/Symbol \242} = ',sprintf("%0.5f",A20),' ", \',"\n", \
      '     ',A21,'*x title " {/Times-Italic A}{/Symbol \242} = ',sprintf("%0.5f",A21),' "'
unset print
update 'Pdep_Cc_lmax100_dl0.001_Op2.dat' 'Pdep_Cc_lmax100_dl0.001_Op2.dat'
set print 'Pdep_Cc_lmax100_dl0.001_Op2.dat'
print "#filename: Pdep_Cc_lmax100_dl0.001_Op2.dat\n", \
      "#Cc\t","A' using exact nu = 0.671696337558\t","P (bar)\n", \
      "#ind. variable\t","(from 3Dvlt_macro.c data fits via 3Dvlt_K0cFind.c)\t","(via my Greywall/Ahlers fit)\n", \
      "1.20\t",A01,"\t",(A01-2.46303)/-0.0281536,"\n", \
      "1.10\t",A02,"\t",(A02-2.46303)/-0.0281536,"\n", \
      "1.06\t",A03,"\t",(A03-2.46303)/-0.0281536,"\n", \
      "1.05\t",A04,"\t",(A04-2.46303)/-0.0281536,"\n", \
      "1.04\t",A05,"\t",(A05-2.46303)/-0.0281536,"\n", \
      "1.03\t",A06,"\t",(A06-2.46303)/-0.0281536,"\n", \
      "1.02\t",A07,"\t",(A07-2.46303)/-0.0281536,"\n", \
      "1.01\t",A08,"\t",(A08-2.46303)/-0.0281536,"\n", \
      "1.00\t",A09,"\t",(A09-2.46303)/-0.0281536,"\n", \
      "0.99\t",A10,"\t",(A10-2.46303)/-0.0281536,"\n", \
      "0.98\t",A11,"\t",(A11-2.46303)/-0.0281536,"\n", \
      "0.97\t",A12,"\t",(A12-2.46303)/-0.0281536,"\n", \
      "0.90\t",A13,"\t",(A13-2.46303)/-0.0281536,"\n", \
      "0.80\t",A14,"\t",(A14-2.46303)/-0.0281536,"\n", \
      "0.70\t",A15,"\t",(A15-2.46303)/-0.0281536,"\n", \
      "0.60\t",A16,"\t",(A16-2.46303)/-0.0281536,"\n", \
      "0.55\t",A17,"\t",(A17-2.46303)/-0.0281536,"\n", \
      "0.50\t",A18,"\t",(A18-2.46303)/-0.0281536,"\n", \
      "0.40\t",A19,"\t",(A19-2.46303)/-0.0281536,"\n", \
      "0.30\t",A20,"\t",(A20-2.46303)/-0.0281536,"\n", \
      "0.20\t",A21,"\t",(A21-2.46303)/-0.0281536,"\n"
unset print
