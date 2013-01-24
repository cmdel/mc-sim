function bsc = BSCall(S0,K,r,sig,T);
%PURPOSE:Returns the European call option price using Black-Scholes
 %-------------------------------------------------------------------------  
 %USAGE : bsc = BSCall(S0,K,r,sig,T)
 %-------------------------------------------------------------------------
 %INPUT :S0 - scalar or vector, price of underlying at time t=0
 %       K - scalar or vector, strike price
 %       r - scalar or vector, continuously compound risk free rate expressed as a
 %       positive decimal number.
 %       sig- scalar or vector, volatility of the underlying(same time units as for r)  
 %       T - scalar or vector, time to maturity (same time units as for r)
 %-------------------------------------------------------------------------
 %RETURN : bsc - scalar or vector, call option price
 %-------------------------------------------------------------------------
echo off;
stau=sqrt(T);
d1=(log(S0./(K.*exp(-r.*T))))./(sig.*stau)+0.5*sig.*stau;
d2=d1-sig.*stau;
bsc=S0.*normcdf(d1)-K.*exp(-r.*T).*normcdf(d2);
echo on;
