
var b c ec lb y ;
varexo u ;

parameters RHO, BETA, M, R, SIGMA ;

model;
ec = c(1);
c = y + b - R*b(-1) ;
lb = 1/c - BETA*R/c(+1) ;
lb = 0;
log(y) = RHO*log(y(-1)) + u ;
end;



shocks;
var u; stderr SIGMA;
end;












