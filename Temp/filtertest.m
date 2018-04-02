% filter test

param.k_n = 2;


Xi = [1,2,4,16,25,36];

% pre-averaged increment
w = ones(1,param.k_n);
pdXi = filter([w,-w],1,Xi);                         % g(x) = min(x,1-x)
pdXi(1:2*param.k_n-1) = 0;                                % edge effect
pdXi = pdXi(2:end);

