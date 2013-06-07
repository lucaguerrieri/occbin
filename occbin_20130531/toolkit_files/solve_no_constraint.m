function [zdata oobase_ Mbase_ ] = ...
    solve_no_constraint(modnam,...
    shockssequence,irfshock,nperiods)

global M_ oo_

errlist = [];

% solve model
eval(['dynare ',modnam,' noclearall'])
oobase_ = oo_;
Mbase_ = M_;
setss

nvars = Mbase_.endo_nbr;
ys_ = oobase_.dr.ys;


[hm1,h,hl1,Jbarmat] = get_deriv(Mbase_,ys_);
cof = [hm1,h,hl1];

[decrulea,decruleb]=get_pq(oobase_.dr);
endog_ = M_.endo_names;
exog_ =  M_.exo_names;




nshocks = size(shockssequence,1);
init = zeros(nvars,1);

wishlist = endog_;
nwishes = size(wishlist,1);


zdata = mkdata(nperiods,decrulea,decruleb,endog_,exog_,wishlist,irfshock,shockssequence);

