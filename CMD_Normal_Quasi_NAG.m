function [quasi,iref,ifail] = GenerateQuasiNormalNAG(xmean,sd,n,dim,qgen)
% Generate Quasi-random 
% Normal variates using 
% different generators:
% qgen=1	Sobol
% qgen=2	Sobol(A659)
% qgen=3	Niederreiter
% qgen=4	Faure

n                           = int64(n);
idim                        = int64(dim);
genid                       = int64(qgen);
xmeanv(1:dim)=xmean;
sdv(1:dim)=sd;
% Skip the first few
% variates in the sequence
iskip                       = int64(1000);

% Initialise the Quasi
% random generator
[iref, ifail]               = g05yl(genid, idim, iskip);


% Generate n quasi random
% normal values
[quasi, iref, ifail]        = g05yj(xmeanv, sdv, n, iref);
