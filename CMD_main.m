%clc;
clear all;
close all;
So=100;
rho=-0.2;
V0=0.04;
xi=0.2;
theta=0.04;
kappa=1.5;
K=100;
T=5.0;
NoSteps=200;  % This is to approximate the trading days in a year for the maturity
NoPaths=50;
lambda=0.0;
r=0.00;
q=0.00; % Non divident stock
NAG=1;     % Set this to the following values to use the NAG libraries 
		   % for Uniform and Normal quasi-random variates
		   % NAG = -1 Disable NAG, use MATLAB's RV generator
		   % NAG = 0  Antithetic variance reduction
		   % NAG = 1  Use NAG libraries for scrambled Quasi-random variates
		   % NAG = 2  USe NAG libraries for unscrambled
% Purturbations of the underlying's price
purt=10;
S=linspace(So-0.8*So,So+0.8*So,purt);
%S=[0.99 1 1.01]
% Start the clock
tic
target = 10;
% Purturb the underlying price
for p = 1:purt
    [P(p), mc(p), err(p,:), V, Si] = CMD_heston_MC(S(p),rho,V0,xi,theta,kappa,K,T,NoSteps,NoPaths,lambda,r,q,NAG);
    if p==target% which prices to show on the plot
    	Vatm = V;
    	Satm = Si;
    	errors = err(:, :);
	end
    HENAG(p) =  s30na('C', K, S(p), T,  xi , kappa, rho, V0, theta, lambda, r, q); 
    deltaHENAG(p) = P(p) -HENAG(p);
end
toc
%% Produce some metrics
fprintf('The standard error is: %g\n',mean(err));
% Create the confidence intervals
alpha = 0.01; % 100-alpha := confidence interbal percientile
Savg = mean(Satm);
errors = [zeros(size(errors,1),1) errors];
err_upper = Savg + norminv(alpha,0,1)*errors(target,:);
err_lower = Savg - norminv(alpha,0,1)*errors(target,:);
% Create random number for the name of the image to have multiple figures
t=floor(rand(1)*100);
figure(t)
time=linspace(0,T,NoSteps+1);
gran=NoPaths/5;
gran=1;
subplot 221;
hold on
plot(time,Satm(1:gran:end,:));
ylabel('Stock price');
xlabel('Time'); 
set(gcf, 'Position', get(0,'Screensize')) % Maximise screen
ylabel('Call Option Price ($)');
plot(time,err_upper, '-.b', 'LineWidth',2);
plot(time,err_lower, '-.b', 'LineWidth',2); 
hold off

subplot 222;
plot(time,Vatm(1:gran:end,:));
ylabel('Variance');
xlabel('Time');

subplot 224;
plot(S,err);
ylabel('Standard error');
xlabel('Underlying price S ($)');
%subplot 225;
%bar(S,deltaHENAG,'r');
%ylabel('Delta to Heston ( MC - NAG )');
%xlabel('Underlying price S ($)');

subplot 223;
plot(S,P, 'Color', 'r');
hold on
plot(S,HENAG);
hold off
ylabel('Call Price');
xlabel('Underlying price S ($)');
if(NAG==-1)
	title(gca,['Call Option prices for the Heston models and ',num2str(T), 'year(s) maturity and MATLAB RV generators']);
elseif(NAG==0)
	title(gca,['Call Option prices for the Heston models and ',num2str(T), 'year(s) maturity and antithetic variance reduction']);
elseif(NAG==1)
	title(gca,['Call Option prices for the Heston models and ',num2str(T), 'year(s) maturity and scattered NAG quasi random variates']);
else
	title(gca,['Call Option prices for the Heston models and ',num2str(T), 'year(s) maturity and NAG quasi random variates']);
end
legend('Monte Carlo Heston','NAG analytical Heston','Location','NorthWest')
 
