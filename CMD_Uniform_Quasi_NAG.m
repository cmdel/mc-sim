function [quasi, irefOut, ifail] = GenerateUniformQuasiNAG(n, genid)
% Generate n Quasi-random Uniform variates vector using the NAG 
% libraries using the following sequence generators
% genid=1	Sobol
% genid=2	Sobol(A659)
% genid=3	Niederreiter
% genid=4	Faure

genid                            = int64(genid);
idim                             = int64(1);
iskip                            = int64(1000);
rcord                            = int64(1);

[iref,  ifail]                   = g05yl(genid, idim, iskip);

[quasi, iref, ifail]             = g05ym(int64(n), rcord, iref);
