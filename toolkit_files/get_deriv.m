function [hm1,h,hl1,j,resid] = get_deriv(M_,ys_)

iy_ = M_.lead_lag_incidence;
it_ = 1;

x = zeros(1,M_.exo_nbr);


% find non-zero columns of hm1
lag_cols = find(iy_(1,:)~=0);
% find non-zero columns of h
con_cols = find(iy_(2,:));
% find non-zero columns of hl1
lea_cols = find(iy_(3,:));


% find number of entries for y vector
ny = length(find(iy_~=0));

% build steady state y
y = ys_(lag_cols);
y = [y;ys_(con_cols)];
y = [y;ys_(lea_cols)];


if ismac
eval(['[resid,g1]=',M_.fname,'_dynamic(y,x, M_.params, it_);']);
else
eval(['[resid,g1]=',M_.fname,'_dynamic(y,x, M_.params, ys_, it_);']);
end


hm1=zeros(M_.endo_nbr);
h = hm1;
hl1 = hm1;
j = zeros(M_.endo_nbr,M_.exo_nbr);


% build hm1
nlag_cols = length(lag_cols);
for i=1:nlag_cols
    hm1(:,lag_cols(i)) = g1(:,i);
end

% build h
ncon_cols = length(con_cols);
for i=1:ncon_cols
    h(:,con_cols(i)) = g1(:,i+nlag_cols);
end

% build hl1
nlea_cols = length(lea_cols);
for i=1:nlea_cols
    hl1(:,lea_cols(i)) = g1(:,i+nlag_cols+ncon_cols);
end


for i = 1:M_.exo_nbr;
    j(:,i) =g1(:,i+ny);
end

