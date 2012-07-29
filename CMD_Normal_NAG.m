function [x,iseedOut,ifail] = GenerateNormalNAG(mu,var, n ,gen)
n = int64(n);
igen = int64(gen);
[igenOut, iseed] = g05kc(igen);
[igen, iseed] = g05kb(igen,iseed);
[x,iseedOut,ifail] = g05la(mu,var,n,igen,iseed);
