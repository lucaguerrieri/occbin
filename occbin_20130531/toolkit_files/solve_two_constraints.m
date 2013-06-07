function [ zdatalinear zdatapiecewise zdatass oo00_  M00_ ] = ...
  solve_two_constraints(modnam_00,modnam_10,modnam_01,modnam_11,...
    constraint1, constraint2,...
    constraint_relax1, constraint_relax2,...
    shockssequence,irfshock,nperiods,curb_retrench,maxiter,init)

global M_ oo_



% solve model
eval(['dynare ',modnam_00,' noclearall nolog'])
oo00_ = oo_;
M00_ = M_;


for i=1:M00_.endo_nbr
  eval([deblank(M00_.endo_names(i,:)) '_ss = oo00_.dr.ys(i); ']);
end
for i = 1:M00_.param_nbr
  eval([M00_.param_names(i,:),'= M00_.params(i);']);
end



eval(['dynare ',modnam_10,' noclearall'])
oo10_ = oo_;
M10_ = M_;

eval(['dynare ',modnam_01,' noclearall'])
oo01_ = oo_;
M01_ = M_;

eval(['dynare ',modnam_11,' noclearall'])
oo11_ = oo_;
M11_ = M_;


% do some error checking

% check inputs
if ~strcmp(M00_.endo_names,M10_.endo_names)
    error([modnam_00,' and ',modnam_10,' need to have exactly the same endogenous variables and they need to be declared in the same order'])
end

if ~strcmp(M00_.exo_names,M10_.exo_names)
    error([modnam_00,' and ',modnam_10,' need to have exactly the same exogenous variables and they need to be declared in the same order'])
end

if ~strcmp(M00_.param_names,M10_.param_names)
    warning(['The parameter list does not match across the files ',modnam_00,' and ',modnam_10])
end


if ~strcmp(M00_.endo_names,M01_.endo_names)
    error([modnam_00,' and ',modnam_01,' need to have exactly the same endogenous variables and they need to be declared in the same order'])
end

if ~strcmp(M00_.exo_names,M01_.exo_names)
    error([modnam_00,' and ',modnam_01,' need to have exactly the same exogenous variables and they need to be declared in the same order'])
end

if ~strcmp(M00_.param_names,M01_.param_names)
    warning(['The parameter list does not match across the files ',modnam_00,' and ',modnam_01])
end


if ~strcmp(M00_.endo_names,M11_.endo_names)
    error([modnam_00,' and ',modnam_11,' need to have exactly the same endogenous variables and they need to be declared in the same order'])
end

if ~strcmp(M00_.exo_names,M11_.exo_names)
    error([modnam_00,' and ',modnam_11,' need to have exactly the same exogenous variables and they need to be declared in the same order'])
end

if ~strcmp(M00_.param_names,M11_.param_names)
    warning(['The parameter list does not match across the files ',modnam_00,' and ',modnam_11])
end





nvars = M00_.endo_nbr;
zdatass = oo00_.dr.ys;


[hm1,h,hl1,Jbarmat] = get_deriv(M00_,zdatass);
cof = [hm1,h,hl1];


M10_.params = M00_.params;
[hm1,h,hl1,Jbarmat10,resid] = get_deriv(M10_,zdatass);
cof10 = [hm1,h,hl1];
Dbarmat10 = resid;

M01_.params = M00_.params;
[hm1,h,hl1,Jbarmat01,resid] = get_deriv(M01_,zdatass);
cof01 = [hm1,h,hl1];
Dbarmat01 = resid;

M11_.params = M00_.params;
[hm1,h,hl1,Jbarmat11,resid] = get_deriv(M11_,zdatass);
cof11 = [hm1,h,hl1];
Dbarmat11 = resid;



[decrulea,decruleb]=get_pq(oo00_.dr);
endog_ = M00_.endo_names;
exog_ =  M00_.exo_names;


% processes the constrain so as to uppend a suffix to each
% endogenous variables
constraint1_difference = process_constraint(constraint1,'_difference',M00_.endo_names,0);

% when the last argument in process_constraint is set to 1, the
% direction of the inequality in the constraint is inverted
constraint_relax1_difference = process_constraint(constraint_relax1,'_difference',M00_.endo_names,0);


% processes the constrain so as to uppend a suffix to each
% endogenous variables
constraint2_difference = process_constraint(constraint2,'_difference',M00_.endo_names,0);

% when the last argument in process_constraint is set to 1, the
% direction of the inequality in the constraint is inverted
constraint_relax2_difference = process_constraint(constraint_relax2,'_difference',M00_.endo_names,0);



nshocks = size(shockssequence,1);




if ~exist('init')
    init = zeros(nvars,1);
end

if ~exist('maxiter')
    maxiter = 20;
end

if ~exist('curb_retrench')
    curb_retrench = 0;
end

init_orig = init;






zdatapiecewise = zeros(nperiods,nvars);


regime1(1) =1;
regime2(1) =1;
regimestart1 =1;
regimestart2 =1;
violvecbool = zeros(nperiods+1,2);  % This sets the first guess for when
% the constraints are going to hold.
% The variable is a boolean with two
% columns. The first column refers to
% constraint1; the second to
% constraint2.
% Each row is a period in time.
% If the boolean is true it indicates
% the relevant constraint is expected
% to evaluate to true.
% The default initial guess is
% consistent with the base model always
% holding -- equivalent to the linear
% solution.

wishlist = endog_;
nwishes = size(wishlist,1);
for ishock = 1:nshocks
    
    
    changes=1;
    iter = 0;
    
    while (changes & iter<maxiter)
        iter = iter +1;
        
        % analyse violvec and isolate contiguous periods in the other
        % regime.
        [regime1 regimestart1]=map_regime(violvecbool(:,1));
        [regime2 regimestart2]=map_regime(violvecbool(:,2));
        
        
        [zdatalinear]=mkdatap_anticipated_2constraints(nperiods,decrulea,decruleb,...
            cof,Jbarmat,...
            cof10,Jbarmat10,Dbarmat10,...
            cof01,Jbarmat01,Dbarmat01,...
            cof11,Jbarmat11,Dbarmat11,...
            regime1,regimestart1,...
            regime2,regimestart2,...
            violvecbool,endog_,exog_,...
            irfshock,shockssequence(ishock,:),init);
        
        for i=1:nwishes
            eval([deblank(wishlist(i,:)),'_difference=zdatalinear(:,i);']);
        end
        
        
        
        
        newviolvecbool1 = eval(constraint1_difference);
        relaxconstraint1 = eval(constraint_relax1_difference);
        
        newviolvecbool2 = eval(constraint2_difference);
        relaxconstraint2 = eval(constraint_relax2_difference);
        
        
        
        newviolvecbool = [newviolvecbool1;newviolvecbool2];
        relaxconstraint = [relaxconstraint1;relaxconstraint2];
        
        
        
        % check if changes
        if (max(newviolvecbool(:)-violvecbool(:)>0)) | sum(relaxconstraint(find(violvecbool==1))>0)
            changes = 1;
        else
            changes = 0;
        end
        
        if curb_retrench   % apply Gauss-Sidel idea of slowing down the change in the guess
            % for the constraint -- only relax one
            % period at a time starting from the last
            % one when each of the constraints is true.
            retrench = 0*violvecbool(:);
            if ~isempty(find(relaxconstraint1 & violvecbool(:,1)))
                retrenchpos = max(find(relaxconstraint1 & violvecbool(:,1)));
                retrench(retrenchpos) = 1;
            end
            if ~isempty(find(relaxconstraint2 & violvecbool(:,2)))
                retrenchpos = max(find(relaxconstraint2 & violvecbool(:,2)));
                retrench(retrenchpos+nperiods+1) = 1;
            end
            violvecbool = (violvecbool(:) | newviolvecbool(:))-retrench(:);
        else
            violvecbool = (violvecbool(:) | newviolvecbool(:))-(relaxconstraint(:) & violvecbool(:));
        end
        
        violvecbool = reshape(violvecbool,nperiods+1,2);
        
        
        
    end
    if changes ==1
        display('Did not converge -- increase maxiter')
    end
    
    init = zdatalinear(1,:);
    zdatapiecewise(ishock,:)=init;
    init= init';
    
    % update the guess for constraint violations for next period
    % update is consistent with expecting no additional shocks next period
    violvecbool=[violvecbool(2:end,:);zeros(1,2)];
    
end


zdatapiecewise(ishock+1:end,:)=zdatalinear(2:nperiods-ishock+1,:);

zdatalinear = mkdata(nperiods,decrulea,decruleb,endog_,exog_,wishlist,irfshock,shockssequence,init_orig);

