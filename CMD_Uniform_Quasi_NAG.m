function [quasi, irefOut, ifail] = GenerateUniformQuasiNAG(n, genid)

genid              = int64(genid);
idim               = int64(1);
iskip              = int64(1000);
rcord              = int64(1);
[iref,ifail]       = g05yl(genid,idim,iskip);
[quasi,iref,ifail] = g05ym(int64(n),rcord,iref);
