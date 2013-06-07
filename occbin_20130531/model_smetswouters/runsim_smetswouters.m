%% set inputs for solution below
%  The program produces responses to the shocks selected below under
%  irfshock. Paths for all the endogenous variables in the model selected
%  are produced. VariableName_difference holds the piece-wise linear solution for
%  VariableName.  VariableName_uncdifference holds the linear solution for
%  VariableName.

clear
%setpathdynare4

global M_ oo_

% modname below chooses model
% directory. But simple param choices are made from paramfile in current
% directory.
modnam = 'usmodel';
modnamstar = 'usmodelzlb';


% express the occasionally binding constraint
% in linearized form
% one can use any combination of endogenous variables and parameters 
% declared in the the dynare .mod files
constraint = 'r<-conster'; 
constraint_relax ='rnot>-conster';


% Pick innovation for IRFs
irfshock =char('eb');      % label for innovation for IRFs
                             % needs to be an exogenous variable in the
                             % dynare .mod files







  nper=10;
  
  shockssequence = [
    -15*ones(nper,1)/nper
    ];         % scale factor for simulations
  nperiods = size(shockssequence,1)+50;            %length of IRFs

  


tic
% Solve model, generate model IRFs
[zdatalinear zdatapiecewise zdatass oobase_ Mbase_  ] = ...
  solve_one_constraint(modnam,modnamstar,...
  constraint, constraint_relax,...
  shockssequence,irfshock,nperiods);
toc
                          
% unpack the IRFs                          
for i=1:M_.endo_nbr
  eval([deblank(M_.endo_names(i,:)),'_uncdifference=zdatalinear(:,i);']);
  eval([deblank(M_.endo_names(i,:)),'_difference=zdatapiecewise(:,i);']);
  eval([deblank(M_.endo_names(i,:)),'_ss=zdatass(i);']);
end


% get parameter values
for i=1:Mbase_.param_nbr
    eval([Mbase_.param_names(i,:),'=Mbase_.params(i);'])
end

%% Modify to plot IRFs and decision rules

titlelist = char('r, Policy Interest Rate (percent at AR)','c, Consumption','inve, Investment','y, Output','b, preference shock');
percent = 'PPt dev. from Steady State';
level = 'Level';
ylabels = char(level,percent,percent,percent,percent);

figtitle = '';
line1=[4*(r_difference+conster),c_difference,inve_difference,y_difference,b_difference];
line2=[4*(r_uncdifference+conster),c_uncdifference,inve_uncdifference,y_uncdifference,b_uncdifference];

% Figure 1
legendlist = cellstr(char('Piecewise Linear','Linear (ignores ZLB)'));
figlabel = '';
makechart(titlelist,legendlist,figlabel,ylabels,line1,line2)

