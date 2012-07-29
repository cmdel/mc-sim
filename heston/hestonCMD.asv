%------------------------------------------------------------------%
function [call]=heston(T,S0,K,v0,theta,kappa,sigma,rho,lambda,r)
%------------------------------------------------------------------%
%                                                                  %
% HESTON    Calculates the heston call option price by using       %
%           20-point Gauss-Legendre(GL) Integration in every       %
%           small interval. Values for the GL weights and          %
%           abscissas were taken from:                                       %
% http://www.efunda.com/math/num_integration/findgausslegendre.cfm %
%                                                                  %
%                                                                  %
% The algorithm is based on the closed-form solution given by      %
% Heston in "A closed-form solution for options with stochastic    %
% volatility with applications to bond and currency options (1993)"%          %
%------------------------------------------------------------------%
% USAGE     [call_prices, std_errs] = Heston(S0, r, V0, eta, theta,%
%           kappa, strike, T, M, N)                                %
%                                                                  %
%------------------------------------------------------------------%
% ACKNOWLEDGMENTS:                                                 %
%                                                                  %
%                                                                  %
%------------------------------------------------------------------%
% INPUTS    S0      - Current price of the underlying asset.
%           r       - Annualized continuously compounded risk-free 
%                     rate of return over the life of the option, 
%                     expressed as a positive decimal number.
%           V0      - Current variance of the underlying asset
%           xi      - volatility of volatility
%           theta   - long-term mean
%           kappa   - rate of mean-reversion
%           K       - Vector of strike prices of the option
%           T       - Time to expiration of the option, 
%                     expressed in years
%                                                                  %
%------------------------------------------------------------------%
% OUTPUTS   [Call]  - Prices (i.e., value) of a vector of          %
%                     European call options.                       %
%                                                                  %
%------------------------------------------------------------------%
% EXAMPLE  heston(7,100,100,0.04,0.04,2,0.2,-0.5,20,0.05)
%                                                                  %
%   S0 = 100;
%   r = 0.02;
%   V0 = 0.04;
%   eta = 0.7;
%   theta = 0.06;
%   kappa = 1.5;
%   strike = 85:5:115;
%   T = 0.25;
%                                                                  %
%   M = 2000;  % Number of paths.
%   N = 250;   % Number of time steps per path
%                                                                  %
%   [call_prices, std_errs] = Heston(S0, r, V0, eta, theta, kappa,
%       strike, T, M, N)
%                                                                  %
% call_prices =
%                                                                  %
%   15.9804   11.4069    7.2125    3.9295    2.1213    1.2922    0.8625
%                                                                  %
%                                                                  %
% std_errs =
%                                                                  %
%    0.0198    0.0263    0.0329    0.0367    0.0357    0.0315    0.0268
%                                                                  %
%                                                                  %
%------------------------------------------------------------------%
% Christos  Delivorias - OR MSc Student                            %
% The University of Edinburgh                                      %
% August 2012                                                      %
%------------------------------------------------------------------%

%
% T: Maturity in years
% S0: Asset price
% K: Strike price
% v0: Initial variance
% theta: long variance
% kappa: mean-recersion rate to theta
% sigma: vol of vol
% rho: correlation of dW^s, dW^v
% lambda: market risk value; normally calculated from market data; 
%		  set to 0 for these examples
% r: interest rate
%
% std_err         - Standard deviation of the error due to the Monte-Carlo 
%                     simulation:
%                     (std_err = std(sample)/sqrt(length(sample)))
% Take the log of the Spot price
x0=log(S0);
a=kappa*theta;

% First row in the vector are the weights for the GL integration and the second are the abscissas
GL_20=[0.0176140070678 0.0406014298819 0.0626720482976 ...
        0.0832767415506 0.101930119826 0.118194531969 ...
		0.131688638458 0.142096109327 0.149172986482 ...
        0.15275338714 0.15275338714 0.149172986482 ...
		0.142096109327 0.131688638458 0.118194531969 ...
        0.101930119826 0.0832767415506 0.0626720482976      ...
		0.0406014298819 0.0176140070678; 
		-0.993128599185 -0.963971927278 -0.912234428251 ...
        -0.839116971822 -0.74633190646 -0.636053680727  ...
		-0.510867001951 -0.373706088715 -0.227785851142 ...
        -0.0765265211335 0.0765265211335 0.227785851142   ...
		0.373706088715 0.510867001951 0.636053680727 ...
        0.74633190646 0.839116971822 0.912234428251   ...
		0.963971927278 0.993128599185];
cutoff=100;    
intvl=2;   
           
P1_int=zeros(1,cutoff/intvl);
P2_int=zeros(1,cutoff/intvl);
for j=1:cutoff/intvl
    for k=1:20
        phi=intvl*GL_20(2,k)/2+intvl*(2*j-1)/2;
        P1_int(j)=GL_20(1,k)*(real(exp(-i*phi*log(K))*...
			feval('fInt',1,x0,v0,T,phi,theta,kappa,sigma,...
                    rho,lambda,r,i,a)/(i*phi)))+P1_int(j);
        P2_int(j)=GL_20(1,k)*(real(exp(-i*phi*log(K))*...
			feval('fInt',2,x0,v0,T,phi,theta,kappa,sigma,...
                    rho,lambda,r,i,a)/(i*phi)))+P2_int(j);
    end
    P1_int(j)=intvl/2*P1_int(j);
    P2_int(j)=intvl/2*P2_int(j);
end
P1=0.5+sum(P1_int)/pi;
P2=0.5+sum(P2_int)/pi;
call=S0*P1-K*exp(-r*T)*P2;
 
%%%
function y=fInt(m,x,v,t,phi,theta,kappa,sigma,rho,lambda,r,i,a)
u=[0.5 -0.5];
b= [kappa+lambda-rho*sigma kappa+lambda];
d=sqrt((rho*sigma*phi*1i-b(m)).^2-sigma^2*(2*u(m)*phi*1i-phi.^2));
g=(b(m)-rho*sigma*phi*1i+d)/(b(m)-rho*sigma*phi*1i-d);
C=r*phi*1i*t+a*(b(m)-rho*sigma*phi*1i+d)*t/(sigma^2);
D=(b(m)-rho*sigma*phi*1i+d)/(sigma^2)*((1-exp(d*t))/(1-g*exp(d*t)));
y=exp(C+D.*v+1i*phi*x).*((1-g.*exp(d.*t))./(1-g)).^(-2*a/(sigma^2));
