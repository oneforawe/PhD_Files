#filename: DelokFit.gnu
N = 1
file(n) = 'sprintf("DFinit%d.par",n)'
FIT_LIMIT = 1e-5  # Given finNL.par below,  L=a => 1e-5;  L=b => 1e-6;  etc.
#f(x,y) = a+b*y+c*exp(-(a+b*y)/x)*x/y+d*exp(-2*(a+b*y)/x)*y+(e+f*y+g*y**2)*exp(-3*(a+b*y)/x)
f(x,y) = b*(x+y)
fit f(x,y) 'DelokFit.dat' using 1:2:3:(3) via file(N)
#update 'sprintf("DFinit%d.par",N)' 'sprintf("DFfin%da.par",N)'
#set print 'sprintf("DFfinal%da.par",N)'
#print a,"\t",b,"\t",c,"\t",d,"\t",e,"\t",f,"\t",g
#unset print
