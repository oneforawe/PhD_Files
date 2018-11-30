#filename: fit_vlt_ThermStates_lmax100_A_RhosoRhoVsTvNu.gnu
f01(x) = A01*x**0.6716883522798452
f02(x) = A02*x**0.6716883522798452
f03(x) = A03*x**0.6716883522798452
f04(x) = A04*x**0.6716883522798452
f05(x) = A05*x**0.6716883522798452
f06(x) = A06*x**0.6716883522798452
f07(x) = A07*x**0.6716883522798452
f08(x) = A08*x**0.6716883522798452
f09(x) = A09*x**0.6716883522798452
f10(x) = A10*x**0.6716883522798452
f11(x) = A11*x**0.6716883522798452
f12(x) = A12*x**0.6716883522798452
f13(x) = A13*x**0.6716883522798452
f14(x) = A14*x**0.6716883522798452
f15(x) = A15*x**0.6716883522798452
f16(x) = A16*x**0.6716883522798452
f17(x) = A17*x**0.6716883522798452
f18(x) = A18*x**0.6716883522798452
f19(x) = A19*x**0.6716883522798452
f20(x) = A20*x**0.6716883522798452
f21(x) = A21*x**0.6716883522798452
fit [0:3.6e-8] f01(x) 'vlt_ThermStates_lmax100_A_Cc1.20.dat' using 2:8 via A01
fit [0:3.6e-8] f02(x) 'vlt_ThermStates_lmax100_A_Cc1.10.dat' using 2:8 via A02
fit [0:3.6e-8] f03(x) 'vlt_ThermStates_lmax100_A_Cc1.06.dat' using 2:8 via A03
fit [0:3.6e-8] f04(x) 'vlt_ThermStates_lmax100_A_Cc1.05.dat' using 2:8 via A04
fit [0:3.6e-8] f05(x) 'vlt_ThermStates_lmax100_A_Cc1.04.dat' using 2:8 via A05
fit [0:3.6e-8] f06(x) 'vlt_ThermStates_lmax100_A_Cc1.03.dat' using 2:8 via A06
fit [0:3.6e-8] f07(x) 'vlt_ThermStates_lmax100_A_Cc1.02.dat' using 2:8 via A07
fit [0:3.6e-8] f08(x) 'vlt_ThermStates_lmax100_A_Cc1.01.dat' using 2:8 via A08
fit [0:3.6e-8] f09(x) 'vlt_ThermStates_lmax100_A_Cc1.00.dat' using 2:8 via A09
fit [0:3.6e-8] f10(x) 'vlt_ThermStates_lmax100_A_Cc0.99.dat' using 2:8 via A10
fit [0:3.6e-8] f11(x) 'vlt_ThermStates_lmax100_A_Cc0.98.dat' using 2:8 via A11
fit [0:3.6e-8] f12(x) 'vlt_ThermStates_lmax100_A_Cc0.97.dat' using 2:8 via A12
fit [0:3.6e-8] f13(x) 'vlt_ThermStates_lmax100_A_Cc0.90.dat' using 2:8 via A13
fit [0:3.6e-8] f14(x) 'vlt_ThermStates_lmax100_A_Cc0.80.dat' using 2:8 via A14
fit [0:3.6e-8] f15(x) 'vlt_ThermStates_lmax100_A_Cc0.70.dat' using 2:8 via A15
fit [0:3.6e-8] f16(x) 'vlt_ThermStates_lmax100_A_Cc0.60.dat' using 2:8 via A16
fit [0:3.6e-8] f17(x) 'vlt_ThermStates_lmax100_A_Cc0.55.dat' using 2:8 via A17
fit [0:3.6e-8] f18(x) 'vlt_ThermStates_lmax100_A_Cc0.50.dat' using 2:8 via A18
fit [0:3.6e-8] f19(x) 'vlt_ThermStates_lmax100_A_Cc0.40.dat' using 2:8 via A19
fit [0:3.6e-8] f20(x) 'vlt_ThermStates_lmax100_A_Cc0.30.dat' using 2:8 via A20
fit [0:3.6e-8] f21(x) 'vlt_ThermStates_lmax100_A_Cc0.20.dat' using 2:8 via A21
update 'fit_vlt_ThermStates_lmax100_A_RhosoRhoVsTvNu.par' 'fit_vlt_ThermStates_lmax100_A_RhosoRhoVsTvNu.par'
set print 'fit_vlt_ThermStates_lmax100_A_RhosoRhoVsTvNu.par'
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
      "Put the following in plot_vlt_ThermStates_lmax100_A_RhosoRhoVsTvNu.gnu\n\n", \
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
update 'Pdep_Cc.dat' 'Pdep_Cc.dat'
set print 'Pdep_Cc.dat'
print "#filename: Pdep_Cc.dat\n", \
      "#Cc\t","A' using exact nu = 0.671688352279845\t","P (bar)\n", \
      "#ind. variable\t","(from vlt_ThermStates.c data fits via vlt_K0cFind.c)\t","(via my Greywall/Ahlers fit)\n", \
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
