function [result,ifail] = GenerateUniformNAG(n,gen)
genid = int64(gen);
subid = int64(0);
[state, ifail] = g05kg(genid,subid); % genid 3 for the Mersenne Twister
[state, result, ifail] = g05sa(int64(n), state);
