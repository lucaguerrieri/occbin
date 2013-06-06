nendog = size(M_.endo_names,1);

for i=1:nendog
  eval([ deblank(M_.endo_names(i,:)) '_ss = oo_.dr.ys(i); ']);
end

nparams = size(M_.param_names);

for i = 1:nparams
  eval([M_.param_names(i,:),'= M_.params(i);']);
end

%if numel(oo_.mean)==numel(oo_.dr.ys)
%  for i=1:nendog
%    evalc([ deblank(M_.endo_names(i,:)) '_sss = oo_.mean(i) ' ]);
%  end
%   for i=1:nendog
%     disp([ (M_.endo_names(i,:)) '   ' num2str(oo_.mean(i),6) '   '   num2str(oo_.dr.ys(i),6) ]);
%   end
%end