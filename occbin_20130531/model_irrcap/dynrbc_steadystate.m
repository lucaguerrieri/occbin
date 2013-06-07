
function [ys,check]=dynrbc_steadystate(junk,ys);

global M_

paramfile_dynrbc

nparams = size(M_.param_names,1);
for icount = 1:nparams
    eval(['M_.params(icount) = ',M_.param_names(icount,:),';'])
end


check=0;

kss = ((1/BETA-1+DELTAK)/ALPHA)^(1/(ALPHA-1));
css = -DELTAK*kss +kss^ALPHA;
iss = DELTAK*kss;

k = log(kss);
c = log(css);
i = log(iss);
kprev = k;
lambdak = 0;
a=0;

ys = [ a
c  
i
k
kprev
lambdak ] ;



