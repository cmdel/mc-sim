clc;
close all;
clear;

So=100;
rho=-0.2;
V0=0.04;
xi=0.2;
theta=0.04;
kappa=1.5;
S=100;
K=100;
T=5.0;
NoSteps=T*250;  % This is to approximate the trading days in a year for the maturity
NoPaths=8000;
lambda=0.0;
r=0.05;
q=0.00; % Non divident stock
NAG=-1;     % Set this to the following values to use the NAG libraries 
		   % for Uniform and Normal quasi-random variates
		   % NAG = -1 Disable NAG, use MATLAB's RV generator
		   % NAG = 0  Antithetic variance reduction
		   % NAG = 1  Use NAG libraries for scrambled Quasi-random variates
		   % NAG = 2  USe NAG libraries for unscrambled

tic;
[p,c,err] = CMD_heston_MC(S,rho,V0,xi,theta,kappa,K,T,NoSteps,NoPaths,lambda,r,q,NAG);
t=toc;

fprintf('Time: %g. Error: %g\n', t,err);
h = msgbox('Finished','Finished', 'help');
