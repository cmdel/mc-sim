% Test of different implementations
%
%   dS = \mu S dt + \sqrt(v_t) S dW_t'
%   dv  = kappa (theta - v_t) dt + xi \sqrt(v_t) dW_t''
%
echo off;
addpath('Heston')
clc
clear all
close all
format short

csv = csvread('data.csv');
size = size(csv);
prices = nan(size(2),7);
disp('-------------------------------------------------------------------------------');
for line=1:size(1)
    %Price
    S=csv(line,1);
    K=csv(line,2);
    r=csv(line,3);
    T=csv(line,4);
    V0=csv(line,5);
    %Stochastic Var
    theta=csv(line,6);
    kappa=csv(line,7);
    xi=csv(line,8);
    rho=csv(line,9);
    lambda=csv(line,10);
    M=1000; N=csv(line,11);  % monte carlo settings

    HC = HestonCall(S,K,r,xi,T,V0,kappa,theta,lambda,rho);
    HCQ = HestonCallQuad(kappa,theta,xi,rho,V0,r,T,S,K);
    Christos = hestonChristos(T,S,K,V0,theta,kappa,xi,rho,lambda,r);
    FFT = HestonFFTVanilla(1,S,K,T,r,r,kappa,theta,xi,rho,V0); 
    CMC = Heston(S, r, V0, xi, theta, kappa, K, T, M, N);
    BS = BSCall(S,K,r,V0,T);
    prices(line,:) = [HC,HCQ,Christos,FFT,CMC,BS, 0;];

end
prices(1,end)=10.0932;
textHeader = {'HC','HCQ','Christos','FFT', 'CMC','B-S', 'MaxMC'};


minimums= nan(size(1),1);
deltas = prices;
for line=1:size(1)
    minimums(line,1)=min(prices(line,:));
    deltas(line,:) = prices(line,:)-prices(line,5);
end

disp(textHeader);
disp('-------------------------------------------------------------------------------');
for line=1:size(1)
    fprintf('\nScenario %g\n\n\n',line);
    disp(prices(line,:));
    disp('');
    fprintf('d(%g) ',deltas(line,:));
    fprintf('\n');
    disp('-------------------------------------------------------------------------------');
end