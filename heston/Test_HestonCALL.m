function Test_HestonCALL()
%--------------------------------------------------------------------------
%PURPOSE: tests the formula of Heston's call, plots the difference between
%the BS option price and Heston's option price for different underlying
%prices
%--------------------------------------------------------------------------
%USAGE: Test_HestonCALL()
%--------------------------------------------------------------------------

clc; clear;

K=100;
r=0.04;
T=0.5;
sig=0.01;
rho=0.5;
vt=0.01;
kap=2;
th=0.02;
St=100;
lda=0;

Stv=(70:1:130)';
NSt=size(Stv,1);
Hv=zeros(NSt,1);
BSv=zeros(NSt,1);

for i=1:NSt
    Hv(i)=HestonCall(Stv(i),K,r,sig,T,vt,kap,th,lda,rho);
    BSv(i)=BSCall(Stv(i),K,r,sqrt(th),T);
end
plot(Stv,(Hv-BSv));

