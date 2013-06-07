%% set inputs for solution below
%  The program produces responses to the shocks selected below under
%  irfshock. Paths for all the endogenous variables in the model selected
%  are produced. VariableName_piecewise holds the piece-wise linear solution for
%  VariableName.  VariableName_linear holds the linear solution for
%  VariableName.

global M_ oo_

% modname below chooses model
% directory. But simple param choices are made from paramfile in current
% directory.
mod00 = 'dynrbc';
mod10 = 'dynrbcineg';
mod01 = 'dynrbcirr';
mod11 = 'dynrbcirrineg';

% Model parameter values (when saved in PARAM_EXTRA, they overwrite the 
%  parameter values in paramfile_dynrbc
PSI=0; PSINEG=5; PHI=0.975;
save PARAM_EXTRA PSI PSINEG PHI

% express the occasionally binding constraint
% in linearized form
% one can use any combination of endogenous variables and parameters
% declared in the the dynare .mod files
% constraint1 defines the first constraint
% if constraint1 is true, solution switches to model2
% but if constraint_relax is true, solution reverts to model1

% The -0.0001 value guarantees convergence of solution
% Without that, there would be too much back and forth switching until convergence,
% and convergence would not be achieved
constraint1 = 'i<-0.0001';
constraint_relax1 = 'i>-0.0001';


constraint2 = 'i<PHI-1';
constraint_relax2 ='lambdak<0';


% Pick innovation for IRFs
irfshock =char('erra');      % label for innovation for IRFs
                             % needs to be an exogenous variable in the
                             % dynare .mod files

maxiter = 200;     
tol0 = 1e-8;


% Option=1: impulse responses
% Option=2: random simulation

option=1;

%%%%%%%%%%%%%%%% Inputs stop here %%%%%%%%%%%%%%%%%%%%%
%% 

if option==1
  nper=1;
  
  shockssequence = [
    zeros(9,1)-0.0001
    -0.01*ones(nper,1)/nper
    zeros(39,1)
    -0.02*ones(nper,1)/nper
    zeros(39,1)
    0.01*ones(nper,1)/nper
    zeros(39,1)
    0.02*ones(nper,1)/nper
    zeros(39,1)
    ];        
  nperiods = size(shockssequence,1);           
  
end

if option==2
   nperiods = 100;
   randn('seed',1);
   shockssequence = 1*randn(nperiods,1)*0.02 ; 
end



% Solve model, generate model IRFs

[zdatalinear zdatapiecewise zdatass oobase_ Mbase_ ] = ...
          solve_two_constraints(mod00,mod10,mod01,mod11,...
                               constraint1, constraint2,...
                               constraint_relax1, constraint_relax2,...
                               shockssequence,irfshock,nperiods);
                          

                          
% unpack the IRFs                          
for i=1:Mbase_.endo_nbr
  eval([deblank(Mbase_.endo_names(i,:)),'_linear=zdatalinear(:,i);']);
  eval([deblank(Mbase_.endo_names(i,:)),'_piecewise=zdatapiecewise(:,i);']);
  eval([deblank(Mbase_.endo_names(i,:)),'_ss=zdatass(i);']);
end


%% Modify to plot IRFs and decision rules


titlelist = char('c','lambdak','k','i','a');
percent = 'Percent';
value = 'value';
ylabels = char(percent,value,percent,percent,percent);
figtitle = 'Simulated variables';
legendlist = cellstr(char('Piecewise Linear','Linear'));

line1=100*[c_piecewise,lambdak_ss/100+lambdak_piecewise/100,(k_piecewise),i_piecewise,a_piecewise];
line2=100*[c_linear,lambdak_ss/100+lambdak_linear/100,(k_linear),i_linear,a_linear ];

makechart9(titlelist,legendlist,figtitle,-1000,ylabels,line1,line2);


