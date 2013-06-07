%% set inputs for solution below
%  The program produces responses to the shocks selected below under
%  irfshock. Paths for all the endogenous variables in the model selected
%  are produced. VariableName_difference holds the piece-wise linear solution for
%  VariableName.  VariableName_uncdifference holds the linear solution for
%  VariableName.

clear 
global M_ oo_

modnam = 'dynrbc';
modnamstar = 'dynrbcirr_i';

% Constraint (see notes 1 and 2 in readme.pdf file)
constraint = 'i<log(PHII)';
constraint_relax ='lambdak<0';

% Pick innovation for IRFs
irfshock =char('erra');      % label for innovation for IRFs

shockssequence = zeros(60,1);
shockssequence(10) = 0.02; 
shockssequence(30) = -0.02; 

maxiter = 10;

nperiods=100;
    
[zdatalinear zdatapiecewise zdatass oobase_ Mbase_  ] = ...
  solve_one_constraint(modnam,modnamstar,...
  constraint, constraint_relax,...
  shockssequence,irfshock,nperiods,maxiter);
    
% unpack the IRFs
for i=1:Mbase_.endo_nbr
  eval([deblank(Mbase_.endo_names(i,:)),'_l=zdatalinear(:,i);']);
  eval([deblank(Mbase_.endo_names(i,:)),'_p=zdatapiecewise(:,i);']);
  eval([deblank(Mbase_.endo_names(i,:)),'_ss=zdatass(i);']);
end


% define inputs for plotting
titlelist = char('c (consumption)','k (capital)','i (investment)','a (tfp)','lambdak (multiplier)');
percent = 'Percent';
level = 'Level';
ylabels = char(percent,percent,percent,percent,level);
figtitle = '';
legendlist = cellstr(char('Piecewise Linear','Linear'));
line1=100*[c_p,k_p,i_p,a_p,(lambdak_p+lambdak_ss)/100];
line2=100*[c_l,k_l,i_l,a_l,(lambdak_l+lambdak_ss)/100];

% create plots
makechart(titlelist,legendlist,figtitle,ylabels,line1,line2)
