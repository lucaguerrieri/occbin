function [zdatalinear zdatapiecewise zdatass oobase_ Mbase_  ] = ...
    solve_one_constraint(modnam,modnamstar,...
    constraint, constraint_relax,...
    shockssequence,irfshock,nperiods,maxiter,init)

global M_ oo_

errlist = [];

% solve model
eval(['dynare ',modnam,' noclearall nolog '])
oobase_ = oo_;
Mbase_ = M_;
setss

eval(['dynare ',modnamstar,' noclearall nolog '])
oostar_ = oo_;
Mstar_ = M_;


% check inputs
if ~strcmp(Mbase_.endo_names,Mstar_.endo_names)
    error('The two .mod files need to have exactly the same endogenous variables declared in the same order')
end

if ~strcmp(Mbase_.exo_names,Mstar_.exo_names)
    error('The two .mod files need to have exactly the same exogenous variables declared in the same order')
end

if ~strcmp(Mbase_.param_names,Mstar_.param_names)
    warning('The parameter list does not match across .mod files')
end

% ensure that the two models have the same parameters
% use the parameters for the base model.
Mstar_.params = Mbase_.params;

nvars = Mbase_.endo_nbr;
zdatass = oobase_.dr.ys;


[hm1,h,hl1,Jbarmat] = get_deriv(Mbase_,zdatass);
cof = [hm1,h,hl1];

[hm1,h,hl1,Jstarbarmat,resid] = get_deriv(Mstar_,zdatass);
cofstar = [hm1,h,hl1];
Dstarbarmat = resid;

[decrulea,decruleb]=get_pq(oobase_.dr);
endog_ = M_.endo_names;
exog_ =  M_.exo_names;


% processes the constrain so as to uppend a suffix to each
% endogenous variables
constraint_difference = process_constraint(constraint,'_difference',Mbase_.endo_names,0);

constraint_relax_difference = process_constraint(constraint_relax,'_difference',Mbase_.endo_names,0);



nshocks = size(shockssequence,1);


if ~exist('init')
    init = zeros(nvars,1);
end

if ~exist('maxiter')
    maxiter = 20;
end

if ~exist('nperiods')
    nperiods = 100;
end

init_orig = init;


zdatapiecewise = zeros(nperiods,nvars);
wishlist = endog_;
nwishes = size(wishlist,1);

violvecbool = zeros(nperiods+1,1);
for ishock = 1:nshocks
    
    changes=1;
    iter = 0;
    
    
    while (changes & iter<maxiter)
        iter = iter +1;
        
        
        [regime regimestart]=map_regime(violvecbool);
        
        
        [zdatalinear]=mkdatap_anticipated(nperiods,decrulea,decruleb,...
            cof,Jbarmat,cofstar,Jstarbarmat,Dstarbarmat,...
            regime,regimestart,violvecbool,...
            endog_,exog_,irfshock,shockssequence(ishock,:),init);
        
        for i=1:nwishes
            eval([deblank(wishlist(i,:)),'_difference=zdatalinear(:,i);']);
        end
        
        
        
        newviolvecbool = eval(constraint_difference);
        relaxconstraint = eval(constraint_relax_difference);
        
        
        
        % check if changes
        if (max(newviolvecbool-violvecbool>0)) | sum(relaxconstraint(find(violvecbool==1))>0)
            changes = 1;
        else
            changes = 0;
        end
        
        
        violvecbool = (violvecbool|newviolvecbool)-(relaxconstraint & violvecbool);
        
        
    end
    
    init = zdatalinear(1,:);
    zdatapiecewise(ishock,:)=init;
    init= init';
    
    % reset violvecbool for next period -- consistent with expecting no
    % additional shocks
    violvecbool=[violvecbool(2:end);0];
    
end


zdatapiecewise(ishock+1:end,:)=zdatalinear(2:nperiods-ishock+1,:);


zdatalinear = mkdata(max(nperiods,size(shockssequence,1)),decrulea,decruleb,endog_,exog_,wishlist,irfshock,shockssequence,init_orig);

if changes ==1
    display('Did not converge -- increase maxiter')
end
