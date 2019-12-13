#filename: calc_Ferer_XvsP_2019.gnu
reset
set output 'calc_Ferer_XvsP_2019.dat'

# The following values and parameters were taken from
# plot_TheoryVsBrooksDonnelly_U0VsP.gnu
# plot_TheoryVsGlaberson_a0VsP.gnu
# 3Dvlt_macro_cp.c version 9  (constants same in version 10, OK!!)
# and
# (for B', experimental/exp and theoretical/loop)
# HeatCapFits3.dat
# from file:///home/zero/Documents/School/UCLA/Research/Williams/programs/results/plots_macro_cp/3D_09_Thesis/  (OLD VERSION!! NEED THE NEW ONE!!)
#
# (NEW FOR 2019 calculations)
# NEW:
# HeatCapFits.dat
# from file:///home/zero/Documents/School/UCLA/Research/Williams/programs/results/plots_macro_cp/3D_14_itlogflow_dl0.0001_Thesis/
#



m = 6.64647585456e-27
NA = 6.02214179e23
R = 8.31447215


# P : Pressure in bars
P1 =  0.050
P2 =  1.646
P3 =  7.328
P4 = 15.031
P5 = 18.180
P6 = 22.533
P7 = 25.868

# B'exp or Be : experimental critical amplitudes for cp ( specific heat capacity, in J/(K.mol) ), corresponding to pressures above
# cp = B'*tau^(-alpha) + cp_critical
Be1 = -388.026880407952
Be2 = -382.149454684935
Be3 = -354.951404759459
Be4 = -337.532647770387
Be5 = -322.721607316235
Be6 = -299.57997589479
Be7 = -306.487231048798
Be(y) = (y == P1 ? Be1 : \
        (y == P2 ? Be2 : \
        (y == P3 ? Be3 : \
        (y == P4 ? Be4 : \
        (y == P5 ? Be5 : \
        (y == P6 ? Be6 : \
        (y == P7 ? Be7 : 0) ))))))

# B'loop or Bl : theoretical (loop) critical amplitudes for cp ( specific heat capacity, in J/(K.mol) ), corresponding to pressures above
# cp = B'*tau^(-alpha) + cp_critical
#Bl1 = -1423.74587326875
#Bl2 = -1432.9178727621
#Bl3 = -1426.76910510232
#Bl4 = -1351.93742251261
#Bl5 = -1309.31833745571
#Bl6 = -1245.86988526681
#Bl7 = -1193.89415030069
Bl1 = -1397.31440791919
Bl2 = -1406.29005656644
Bl3 = -1400.57589816968
Bl4 = -1327.60605121375
Bl5 = -1285.62078993663
Bl6 = -1223.36373708572
Bl7 = -1173.01667705812
Bl(y) = (y == P1 ? Bl1 : \
        (y == P2 ? Bl2 : \
        (y == P3 ? Bl3 : \
        (y == P4 ? Bl4 : \
        (y == P5 ? Bl5 : \
        (y == P6 ? Bl6 : \
        (y == P7 ? Bl7 : 0) ))))))


# A' or A : Parameters & eqn for A(P):
# B' : experimental critical amplitudes for cp ( specific heat capacity, in J/(K.mol) ), corresponding to pressures above
# cp = B'*tau^(-alpha) + cp_critical
A0 = 2.463
A1 = -0.02815
A(y) = A0 + A1*y

# Tc : Parameters & eqn for Tc(P):
# Tc expressed in kelvins, P must be expressed in bars
# HEY, shouldn't Tc0 be 2.1768 K, for P=0?
Tc0 = 2.17349425585161
Tc1 = -0.00982499579394534
Tc2 = -0.000118194448444384
Tc3 = -4.36914591522034e-07
Tc4 = 7.39407378262721e-09
Tc(y)  = Tc0 + Tc1*y + Tc2*y**2 + Tc3*y**3 + Tc4*y**4

# K0c : Parameters & eqn for K0c(P):
k00 = 0.295359928991731
k01 = 0.00508413174676733
k02 = 5.15885451496554e-05
k03 = 6.69186316462475e-07
k04 = -2.74622551950506e-08
k05 = 2.92854219547921e-09
k06 = -8.80012434467914e-11
k07 = -1.93460846518787e-13
k08 = 7.13136752136051e-14
k09 = -1.47121339074112e-15
k10 = 9.60439029245856e-18
K0c(y)  = k00 + k01*y + k02*y**2 + k03*y**3 + k04*y**4 + k05*y**5 + k06*y**6 + k07*y**7 + k08*y**8 + k09*y**9 + k10*y**10

# rho : Parameters & eqn for rho(T,P):
# rho expressed in (kg m^-3), T,P must be expressed in kelvins,bars
R00 = 145.145109496329
R10 = -0.0976539693059151
R20 = 0.334163407001684
R30 = -0.446930785976304
R40 = 0.181879478545246
R01 = 1.74477604495583
R11 = -0.0919538993179052
R21 = 0.179844560873926
R31 = -0.133606331352667
R41 = 0.0410225514249919
R02 = -0.0491655379690169
R12 = 0.00710698898070406
R22 = -0.00823054225495917
R32 = 0.000609542602247143
R42 = 0.00114916775392305
R03 = 0.0013415037643754
R13 = -0.000362007479155809
R23 = 0.000358809384119286
R33 = 6.48183954357527e-05
R43 = -0.000104112551302631
R04 = -1.69907294147191e-05
R14 = 5.53820368251513e-06
R24 = -3.15773411117433e-06
R34 = -4.99967306908062e-06
R44 = 3.41331223468399e-06
rho(x,y) = R00      + R10*x      + R20*x**2      + R30*x**3      + R40*x**4      + R01*y    + R11*x*y    + R21*x**2*y    + R31*x**3*y    + R41*x**4*y    + R02*y**2 + R12*x*y**2 + R22*x**2*y**2 + R32*x**3*y**2 + R42*x**4*y**2 + R03*y**3 + R13*x*y**3 + R23*x**2*y**3 + R33*x**3*y**3 + R43*x**4*y**3 + R04*y**4 + R14*x*y**4 + R24*x**2*y**4 + R34*x**3*y**4 + R44*x**4*y**4
rhoc(y)  = rho(Tc(y),y)

Vm(x,y) = NA*m/rho(x,y)
Vc(y)   = Vm(Tc(y),y)

#experimental:
Xe(y) = (1/R) * ( Be(y) / Vc(y) ) * ( (Tc(y)*Vc(y)) / A(y) )**3

#theoretical/loop:
Xl(y) = (1/R) * ( Bl(y) / Vc(y) ) * ( (Tc(y)*Vc(y)) / A(y) )**3


#P	B'	A'	Tc	rhoc	Vc

set print 'calc_Ferer_XvsP_2019.dat'
print "filename: calc_Ferer_XvsP_2019.dat\n\n", \
      "P","\t","B'exp","\t","B'loop","\t","A'","\t","Tc","\t","rhoc","\t","Vc","\t","Xexp","\t","Xloop"
P = P1
print P,"\t",Be(P),"\t",Bl(P),"\t",A(P),"\t",Tc(P),"\t",rhoc(P),"\t",Vc(P),"\t",Xe(P),"\t",Xl(P)
P = P2
print P,"\t",Be(P),"\t",Bl(P),"\t",A(P),"\t",Tc(P),"\t",rhoc(P),"\t",Vc(P),"\t",Xe(P),"\t",Xl(P)
P = P3
print P,"\t",Be(P),"\t",Bl(P),"\t",A(P),"\t",Tc(P),"\t",rhoc(P),"\t",Vc(P),"\t",Xe(P),"\t",Xl(P)
P = P4
print P,"\t",Be(P),"\t",Bl(P),"\t",A(P),"\t",Tc(P),"\t",rhoc(P),"\t",Vc(P),"\t",Xe(P),"\t",Xl(P)
P = P5
print P,"\t",Be(P),"\t",Bl(P),"\t",A(P),"\t",Tc(P),"\t",rhoc(P),"\t",Vc(P),"\t",Xe(P),"\t",Xl(P)
P = P6
print P,"\t",Be(P),"\t",Bl(P),"\t",A(P),"\t",Tc(P),"\t",rhoc(P),"\t",Vc(P),"\t",Xe(P),"\t",Xl(P)
P = P7
print P,"\t",Be(P),"\t",Bl(P),"\t",A(P),"\t",Tc(P),"\t",rhoc(P),"\t",Vc(P),"\t",Xe(P),"\t",Xl(P)


unset print














# Adjusting units to match Ferer:
Vc(y) = (10**6)*Vm(Tc(y),y)

#experimental:
Xe(y) = (10**-4) * (1/R) * ( Be(y) / Vc(y) ) * ( (Tc(y)*Vc(y)) / A(y) )**3

#theoretical/loop:
Xl(y) = (10**-4) * (1/R) * ( Bl(y) / Vc(y) ) * ( (Tc(y)*Vc(y)) / A(y) )**3



set print 'calc_Ferer_XvsP_2019_4Table.dat'
print "filename: calc_Ferer_XvsP_2019_4Table.dat\n", \
      "(derived from the SHIFTED POWER FORM FIT)\n\n", \
      "P\tX1\tX2\tXexp\tXloop\n", \
      "(bar)\t((R/cm3)(Kcm3/mol)3 ...\n"
P = P1
print " $0$&$05$ & 2.06 & 3.11 & $", sprintf("%3.2f",(Xe(P))), "$ & $", sprintf("%3.2f",(Xl(P))), '$  \\'

P = P2
print " $1$&$65$ & 2.02 & 3.07 & $", sprintf("%3.2f",(Xe(P))), "$ & $", sprintf("%3.2f",(Xl(P))), '$  \\'

P = P3
print " $7$&$33$ & 1.99 & 3.01 & $", sprintf("%3.2f",(Xe(P))), "$ & $", sprintf("%3.2f",(Xl(P))), '$  \\'

P = P4
print "$15$&$03$ & 2.02 & 3.10 & $", sprintf("%3.2f",(Xe(P))), "$ & $", sprintf("%3.2f",(Xl(P))), '$  \\'

P = P5
print "$18$&$18$ & 1.99 & 3.10 & $", sprintf("%3.2f",(Xe(P))), "$ & $", sprintf("%3.2f",(Xl(P))), '$  \\'

P = P6
print "$22$&$53$ & 1.96 & 3.12 & $", sprintf("%3.2f",(Xe(P))), "$ & $", sprintf("%3.2f",(Xl(P))), '$  \\'

P = P7
print "$25$&$87$ & 2.19 & 3.47 & $", sprintf("%3.2f",(Xe(P))), "$ & $", sprintf("%3.2f",(Xl(P))), '$  \bstrut\\'

unset print















#reset
#set terminal postscript color eps enhanced
#set output 'calc_Ferer_XvsP.ps'

