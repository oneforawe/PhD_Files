#filename: fit_3Dvlt_macro_Tv_Op12_lmax100_dl0.0001_RhosoRhoVsTvNu.gnu

# update these values (from 3Dvlt_K0cFind.c Output) and values below before running
K0c01 = 0.9346549781218052
K0c02 = 0.79674794984036801
K0c03 = 0.70109094206651523
K0c04 = 0.63025269735391265
K0c05 = 0.57534342964893814
K0c06 = 0.53132345421063487
K0c07 = 0.49510828997129153
K0c08 = 0.46469704703948356
K0c09 = 0.43873100720691904
K0c10 = 0.41625279991888953
K0c11 = 0.39656697289619258
K0c12 = 0.37915524531909928
K0c13 = 0.36362286558404033
K0c14 = 0.3496634698826776
K0c15 = 0.33703537107020987
K0c16 = 0.32554514449601912
K0c17 = 0.31503600689857347
K0c18 = 0.3053794236081705
K0c19 = 0.29646893900678273
K0c20 = 0.28821556877326282
K0c21 = 0.28054430897522864

f01(x) = A01*x**0.671688352279845358
f02(x) = A02*x**0.671688352279845358
f03(x) = A03*x**0.671688352279845358
f04(x) = A04*x**0.671688352279845358
f05(x) = A05*x**0.671688352279845358
f06(x) = A06*x**0.671688352279845358
f07(x) = A07*x**0.671688352279845358
f08(x) = A08*x**0.671688352279845358
f09(x) = A09*x**0.671688352279845358
f10(x) = A10*x**0.671688352279845358
f11(x) = A11*x**0.671688352279845358
f12(x) = A12*x**0.671688352279845358
f13(x) = A13*x**0.671688352279845358
f14(x) = A14*x**0.671688352279845358
f15(x) = A15*x**0.671688352279845358
f16(x) = A16*x**0.671688352279845358
f17(x) = A17*x**0.671688352279845358
f18(x) = A18*x**0.671688352279845358
f19(x) = A19*x**0.671688352279845358
f20(x) = A20*x**0.671688352279845358
f21(x) = A21*x**0.671688352279845358
fit [0:3.6e-8] f01(x) '3Dvlt_macro_Tv_Op12_Cc0.20_lmax100_dl0.0001.dat' using 4:10 via A01
fit [0:3.6e-8] f02(x) '3Dvlt_macro_Tv_Op12_Cc0.25_lmax100_dl0.0001.dat' using 4:10 via A02
fit [0:3.6e-8] f03(x) '3Dvlt_macro_Tv_Op12_Cc0.30_lmax100_dl0.0001.dat' using 4:10 via A03
fit [0:3.6e-8] f04(x) '3Dvlt_macro_Tv_Op12_Cc0.35_lmax100_dl0.0001.dat' using 4:10 via A04
fit [0:3.6e-8] f05(x) '3Dvlt_macro_Tv_Op12_Cc0.40_lmax100_dl0.0001.dat' using 4:10 via A05
fit [0:3.6e-8] f06(x) '3Dvlt_macro_Tv_Op12_Cc0.45_lmax100_dl0.0001.dat' using 4:10 via A06
fit [0:3.6e-8] f07(x) '3Dvlt_macro_Tv_Op12_Cc0.50_lmax100_dl0.0001.dat' using 4:10 via A07
fit [0:3.6e-8] f08(x) '3Dvlt_macro_Tv_Op12_Cc0.55_lmax100_dl0.0001.dat' using 4:10 via A08
fit [0:3.6e-8] f09(x) '3Dvlt_macro_Tv_Op12_Cc0.60_lmax100_dl0.0001.dat' using 4:10 via A09
fit [0:3.6e-8] f10(x) '3Dvlt_macro_Tv_Op12_Cc0.65_lmax100_dl0.0001.dat' using 4:10 via A10
fit [0:3.6e-8] f11(x) '3Dvlt_macro_Tv_Op12_Cc0.70_lmax100_dl0.0001.dat' using 4:10 via A11
fit [0:3.6e-8] f12(x) '3Dvlt_macro_Tv_Op12_Cc0.75_lmax100_dl0.0001.dat' using 4:10 via A12
fit [0:3.6e-8] f13(x) '3Dvlt_macro_Tv_Op12_Cc0.80_lmax100_dl0.0001.dat' using 4:10 via A13
fit [0:3.6e-8] f14(x) '3Dvlt_macro_Tv_Op12_Cc0.85_lmax100_dl0.0001.dat' using 4:10 via A14
fit [0:3.6e-8] f15(x) '3Dvlt_macro_Tv_Op12_Cc0.90_lmax100_dl0.0001.dat' using 4:10 via A15
fit [0:3.6e-8] f16(x) '3Dvlt_macro_Tv_Op12_Cc0.95_lmax100_dl0.0001.dat' using 4:10 via A16
fit [0:3.6e-8] f17(x) '3Dvlt_macro_Tv_Op12_Cc1.00_lmax100_dl0.0001.dat' using 4:10 via A17
fit [0:3.6e-8] f18(x) '3Dvlt_macro_Tv_Op12_Cc1.05_lmax100_dl0.0001.dat' using 4:10 via A18
fit [0:3.6e-8] f19(x) '3Dvlt_macro_Tv_Op12_Cc1.10_lmax100_dl0.0001.dat' using 4:10 via A19
fit [0:3.6e-8] f20(x) '3Dvlt_macro_Tv_Op12_Cc1.15_lmax100_dl0.0001.dat' using 4:10 via A20
fit [0:3.6e-8] f21(x) '3Dvlt_macro_Tv_Op12_Cc1.20_lmax100_dl0.0001.dat' using 4:10 via A21
update 'fit_3Dvlt_macro_Tv_Op12_lmax100_dl0.0001_RhosoRhoVsTvNu.par' 'fit_3Dvlt_macro_Tv_Op12_lmax100_dl0.0001_RhosoRhoVsTvNu.par'
set print 'fit_3Dvlt_macro_Tv_Op12_lmax100_dl0.0001_RhosoRhoVsTvNu.par'
print "#filename: fit_3Dvlt_macro_Tv_Op12_lmax100_dl0.0001_RhosoRhoVsTvNu.par\n", \
      "#source: fit_3Dvlt_macro_Tv_Op12_lmax100_dl0.0001_RhosoRhoVsTvNu.gnu\n\n", \
      "A' ( for  Rhos/Rho = A'*(1-T/Tc)^nu  at different Cc's / K0c's / Pressures )\n", \
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
update 'Pdep_CcK0c.dat' 'Pdep_CcK0c.dat'
set print 'Pdep_CcK0c.dat'
print "#filename: Pdep_CcK0c.dat\n", \
      "#source: fit_3Dvlt_macro_Tv_Op12_lmax100_dl0.0001_RhosoRhoVsTvNu.gnu\n", \
      "#Cc\tK0c\t","A' using exact nu = 0.671688352279845358\t","P (bar)\n", \
      "(#ind. variable)\t(from 3Dvlt_K0cFind.c Output)\t","(from 3Dvlt_K0cFind.c via data fits using 3Dvlt_macro.c)\t","(via my Greywall/Ahlers fit)\n", \
      "0.20\t",sprintf("%0.17f",K0c01),"\t",A01,"\t",(A01-2.46303)/-0.0281536,"\n", \
      "0.25\t",sprintf("%0.17f",K0c02),"\t",A02,"\t",(A02-2.46303)/-0.0281536,"\n", \
      "0.30\t",sprintf("%0.17f",K0c03),"\t",A03,"\t",(A03-2.46303)/-0.0281536,"\n", \
      "0.35\t",sprintf("%0.17f",K0c04),"\t",A04,"\t",(A04-2.46303)/-0.0281536,"\n", \
      "0.40\t",sprintf("%0.17f",K0c05),"\t",A05,"\t",(A05-2.46303)/-0.0281536,"\n", \
      "0.45\t",sprintf("%0.17f",K0c06),"\t",A06,"\t",(A06-2.46303)/-0.0281536,"\n", \
      "0.50\t",sprintf("%0.17f",K0c07),"\t",A07,"\t",(A07-2.46303)/-0.0281536,"\n", \
      "0.55\t",sprintf("%0.17f",K0c08),"\t",A08,"\t",(A08-2.46303)/-0.0281536,"\n", \
      "0.60\t",sprintf("%0.17f",K0c09),"\t",A09,"\t",(A09-2.46303)/-0.0281536,"\n", \
      "0.65\t",sprintf("%0.17f",K0c10),"\t",A10,"\t",(A10-2.46303)/-0.0281536,"\n", \
      "0.70\t",sprintf("%0.17f",K0c11),"\t",A11,"\t",(A11-2.46303)/-0.0281536,"\n", \
      "0.75\t",sprintf("%0.17f",K0c12),"\t",A12,"\t",(A12-2.46303)/-0.0281536,"\n", \
      "0.80\t",sprintf("%0.17f",K0c13),"\t",A13,"\t",(A13-2.46303)/-0.0281536,"\n", \
      "0.85\t",sprintf("%0.17f",K0c14),"\t",A14,"\t",(A14-2.46303)/-0.0281536,"\n", \
      "0.90\t",sprintf("%0.17f",K0c15),"\t",A15,"\t",(A15-2.46303)/-0.0281536,"\n", \
      "0.95\t",sprintf("%0.17f",K0c16),"\t",A16,"\t",(A16-2.46303)/-0.0281536,"\n", \
      "1.00\t",sprintf("%0.17f",K0c17),"\t",A17,"\t",(A17-2.46303)/-0.0281536,"\n", \
      "1.05\t",sprintf("%0.17f",K0c18),"\t",A18,"\t",(A18-2.46303)/-0.0281536,"\n", \
      "1.10\t",sprintf("%0.17f",K0c19),"\t",A19,"\t",(A19-2.46303)/-0.0281536,"\n", \
      "1.15\t",sprintf("%0.17f",K0c20),"\t",A20,"\t",(A20-2.46303)/-0.0281536,"\n", \
      "1.20\t",sprintf("%0.17f",K0c21),"\t",A21,"\t",(A21-2.46303)/-0.0281536,"\n"
unset print
