%---------------------------------------------------------------------
% Simple model to compare government spending multipliers at the ZLB
%---------------------------------------------------------------------

clear
close all
set(0,'DefaultLineLineWidth',2)

%---------------------------------------------------------------------
% option=1: solve model disregarding ZLB (calls function solve_no_constraint)
% option=2: take ZLB into account
%---------------------------------------------------------------------
option=2;

%---------------------------------------------------------------------
% To compute multipliers at ZLB, solve two models
% Model 1 is baseline  that takes us at the ZLB
% Model 2 is baseline model that takes us at the ZLB plus G shock
%---------------------------------------------------------------------


for model=1:2
   
  % Pick color of charts
  if model==1;
    this_color='b';
  else
    this_color='r';
  end
  
  % This line to save additional parameters
  ZZZ=NaN; save PARAM_EXTRA ZZZ;
  
  irfshock = char('eps_c','eps_g'); % Shocks we look at: preference and g-shock
  
  % Model 1 is a model with a zero baseline
  if model==1
    baseline0=[ 0     0      0       0        0       0
              0.00  0.00   0.00    0.00     0.00    0.0 ]';
  end
  
  % Model 2 is a model with a ZLB baseline
  if model==2
    baseline0=[ -0.04  -0.04   -0.04    -0.04      0   0
                 0.00   0.00   0.00     0.00     0.0   0.0  ]';
  end
  
  % In both models, there is a positive G shock in period 6
  scenario1=[ 0     0      0       0        0       0
            0.00  0.00   0.00    0.00     0.00    0.01   ]';
  
  T=30;
  
  if option==1
    
    % Here solve model without ZLB
    % Note I am solving model twice
    % First, I solve it under the sequence of shocks described in baseline0
    % Next,  I solve it under the sequence of shocks described in baseline0+scenario1
    modnam = 'dnk';
    [zdata0 oobase_ Mbase_ ] = solve_no_constraint(modnam,baseline0,irfshock,T) ;
    [zdata1 oobase_ Mbase_ ] = solve_no_constraint(modnam,baseline0+scenario1,irfshock,T) ;
    
  elseif option==2
    
    modnam = 'dnk';
    modnamstar = 'dnk_zlb';
     
    constraint = 'r<-(1/BETA-1)';
    constraint_relax ='rnot>-(1/BETA-1)';
    
    % First time we solve model only with baseline shocks
    % zdata0 contains impulse responses
    [zdatau zdata0 zdatass oobase_ Mbase_] = ...
      solve_one_constraint(modnam,modnamstar,...
      constraint, constraint_relax,...
      baseline0,irfshock,T,20);
    
    % Second time we solve model with baseline shocks and scenario 
    % zdata1 contains impulse responses
    [zdatau zdata1 zdatass oobase_ Mbase_ ] = ...
      solve_one_constraint(modnam,modnamstar,...
      constraint, constraint_relax,...
      baseline0+scenario1,irfshock,T,20);
    
  end
  
  % Note that we compute impulse responses in deviation from baseline0
  % In model=1, baseline0 has a no negative preference shock 
  % In model=2, baseline0 has a negative preference shock that takes economy to ZLB
  for i=1:Mbase_.endo_nbr
    eval([deblank(Mbase_.endo_names(i,:)),'=zdata1(:,i)-zdata0(:,i);']);
  end
  
  
  scale=100;
  subplot(3,2,1)
  plot(scale*y,'color',this_color); hold on
  title('Output')
  ylabel('percent deviation from baseline')

  subplot(3,2,2)
  plot(scale*a_g,'color',this_color); hold on
  title('G/Y')
  ylabel('% of GDP, deviation from baseline')
  xlabel('quarters')
  
  subplot(3,2,3)
  plot(scale*4*r,'color',this_color); hold on
  ylabel('ppoints deviation from baseline, annualized')
  xlabel('quarters')
  title('Interest rate')
  
  subplot(3,2,4)
  plot(scale*c,'color',this_color); hold on
  ylabel('percent deviation from baseline')
  xlabel('quarters')
  title('Consumption')
  
  subplot(3,2,5)
  plot(scale*ik,'color',this_color); hold on
  ylabel('percent deviation from baseline')
  xlabel('quarters')
  title('Investment')
  
  legend('Model with a zero baseline','Model with negative baseline and potential ZLB')
  
  
end


