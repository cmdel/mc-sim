% Test of different implementations
%
%   dS = \mu S dt + \sqrt(v_t) S dW_t'
%   dv  = kappa (theta - v_t) dt + xi \sqrt(v_t) dW_t''
%
clc
clear all
close all
format short
addpath('../')

tic;
points=21;
fprintf('Reading file....\n\n');
csv = csvread('data.csv');
size = size(csv);
scenarios = nan(size(2),1);

for line=1:size(1)
    fprintf('Begin Scenario %g....\n\n',line);
    %Price
    S=linspace(csv(line,1)-0.8*csv(line,1),csv(line,1)+0.8*csv(line,1),points);
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
    q=csv(line,11);
    steps=csv(line,12); 
    paths=csv(line,13);  % monte carlo settings
    for y = 1:points 
        HE_NAG(y) = s30na('C' , K , S(y) , T , xi , kappa, rho, V0, theta, lambda, r, q);
        HE_MC_CMD(y) = CMD_heston_MC(S(y),rho,V0,xi,theta, kappa, K, T, steps, paths, lambda, r, q, 1);
        BS(y) = BSCall(S(y),K,r,V0,T);
        BS_C(y) = s30aa('C' , K , S(y) , T , V0 , r , q);
    end
    MaxMC = [0.0 0.0 0.0 0.003 0.007 0.058 0.291 0.979 2.776 5.734 10.182 15.819 22.211 30.012 37.282 45.113 52.777 60.841 68.507 76.651 83.557];
    figure(line);
    plot(S, HE_NAG,'-s','Color','r')
    title(gca,['Call Option prices for the Heston models and maturity T = ',num2str(T),' year(s)']);
    hold on
    plot(S, BS_C,'-d','Color','k')
    plot(S, BS,'--','Color','c')
    plot(S, HE_MC_CMD,'-d','Color','m')
    if line == 1
        plot(S, MaxMC,'-p','Color','y')
    end 
    ylabel('Call Option Price ($)');
    xlabel('Underlying price S ($)');
    ylim([-10 130]);
    legend('Heston - NAG','Black-Scholes - NAG','Black-Scholes - Analytical','Heston - Quasi-MC with QE','MaxMC-10^6 paths','Location','NorthWest')
    hold off
end

toc;
