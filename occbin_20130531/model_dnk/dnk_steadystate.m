
function [ys,check]=dnk_steadystate(junk,ys)

global M_

paramfile_dnk

nparams = size(M_.param_names,1);
for icount = 1:nparams
    eval(['M_.params(icount) = ',M_.param_names(icount,:),';'])
end


check=0;

% Interest rates
r = 1 / BETA ;
rk = r - (1-DK) ;



K_TO_Y = BETA*ALPHA/(1-BETA*(1-DK))/XP_SS ;
C_TO_Y = (XP_SS-1+(r-1)*K_TO_Y*XP_SS+(1-ALPHA))/XP_SS - GBAR  ;
n = ((1-ALPHA)/(C_TO_Y)/XP_SS/XW_SS/TAU)^(1/(1+ETA)) ;
y = (n) *  K_TO_Y^(ALPHA/(1-ALPHA)) ;
k = K_TO_Y*y ;
c = C_TO_Y*y  ;
uc = 1/c;
un = TAU*n^ETA;
w = (1-ALPHA)*y/XP_SS/n ;

ik = DK * k   ;

dp = 1 ;


zk = 1;


xw = XW_SS;
xp = XP_SS;

mack = 1;
vk = 1;

a_z=1;
a_g=GBAR ;
a_c=1;

b=(1-ETAXG)/(ETAXB+1-r)*a_g ;
tax=ETAXB*b+ETAXG*a_g;
p = 1;

% b, tax and a_g are in levels
xxx = [ ...
a_c
exp(a_g)
a_z
exp(b)
c
dp
ik 
k
mack
n
p
r
rk
r
exp(tax)
uc
un
vk
w
xp
xw
y
zk  ];



logxxx = log(xxx) ;


ys = logxxx ;