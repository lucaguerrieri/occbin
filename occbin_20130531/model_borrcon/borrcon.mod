
var b c ec lb y ;
varexo u ;

parameters RHO, BETA, M, R, SIGMA  ;



model;
ec = c(1);
c = y + b - R*b(-1) ;
b =  M*y ;
lb = 1/c - BETA*R/c(+1) ;
log(y) = RHO*log(y(-1)) + u ;
end;





shocks;
var u; stderr SIGMA;
end;

steady;




stoch_simul(order=1,nocorr,nomoments,irf=0);









