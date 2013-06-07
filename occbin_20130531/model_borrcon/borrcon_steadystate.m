
function [ys,check]=borrcon_steadystate(junk,ys);

global M_

paramfile_borrcon

nparams = size(M_.param_names,1);
for icount = 1:nparams
    eval(['M_.params(icount) = ',M_.param_names(icount,:),';'])
end


check=0;

b=M;
c=1+M-R*M;
ec=c;
lb=(1-BETA*R)/c;
y=1;


ys = [ b
c  
ec
lb
y ] 



