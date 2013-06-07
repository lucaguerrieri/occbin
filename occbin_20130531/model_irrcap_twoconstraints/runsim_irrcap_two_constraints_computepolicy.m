%% set inputs for solution below
%  The program produces responses to the shocks selected below under
%  irfshock. Paths for all the endogenous variables in the model selected
%  are produced. VariableName_piecewise holds the piece-wise linear solution for
%  VariableName.  VariableName_linear holds the linear solution for
%  VariableName.

clear
global M_ oo_

% modname below chooses model
% directory. But simple param choices are made from paramfile in current
% directory.
mod00 = 'dynrbc';
mod10 = 'dynrbcineg';
mod01 = 'dynrbcirr';
mod11 = 'dynrbcirrineg';


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

maxiter = 20;
nperiods = 100;

a_vec = linspace(-0.03,0.03,20);
% k_vec = linspace(-0.05,0.05,3);
k_vec = linspace(0,0,1);

for index_a=1:numel(a_vec)
  for index_k=1:numel(k_vec)
    
    % Order of variables is a, c, i, k, lambdak
    initcon=[ 0 0 0 k_vec(index_k) 0 ]
    shockssequence=[ a_vec(index_a) ]
    
    curb_retrench = 0;
    [zdatalinear zdatapiecewise zdatass oo_ M_ ] = ...
      solve_two_constraints(mod00,mod10,mod01,mod11,...
      constraint1, constraint2,...
      constraint_relax1, constraint_relax2,...
      shockssequence,irfshock,nperiods,curb_retrench,maxiter,initcon);
    
    % unpack the IRFs
    for i=1:M_.endo_nbr
      eval([deblank(M_.endo_names(i,:)),'_linear=zdatalinear(:,i);']);
      eval([deblank(M_.endo_names(i,:)),'_piecewise=zdatapiecewise(:,i);']);
      eval([deblank(M_.endo_names(i,:)),'_ss=zdatass(i);']);
    end
    
    
    kdec(index_a,index_k)=k_piecewise(1)+k_ss;
    idec(index_a,index_k)=i_piecewise(1)+i_ss;
    cdec(index_a,index_k)=c_piecewise(1)+c_ss;
    astate(index_a,index_k)=a_vec(index_a)+a_ss;
    kstate(index_a,index_k)=k_vec(index_k)+a_ss;
    
  end
end

figure
for ik=1:numel(k_vec)
  subplot(2,1,1)
  plot(astate(:,ik),idec(:,ik)-i_ss,'color',[ ik/numel(k_vec) 0 0]); hold on
  xlabel('Technology deviation from ss')  
  ylabel('Log Investment minus ss')  
  title('Policy function for investment as a function of technology')
  subplot(2,1,2)
  plot(astate(:,ik),(idec(:,ik)-i_ss)./astate(:,ik),'color',[ ik/numel(k_vec) 0 0]); hold on
  xlabel('Technology deviation from ss')  
  ylabel('Elasticity of inv. to tec.')
end
if numel(k_vec)==2
legend(['k(t-1)=' num2str(k_vec(1))],['k(t-1)=' num2str(k_vec(2))]) 
elseif numel(k_vec)==3
legend(['k(t-1)=' num2str(k_vec(1))],['k(t-1)=' num2str(k_vec(2))],['k(t-1)=' num2str(k_vec(3))])  
end
