%% Monte Carlo Simulation of the Heston Model
% Calculate the value of a European call vanilla option
% using the Heston model for stochastic volatility

%% USAGE
% [call,std_err]=MCheston(So, r, V0, xi, theta, kappa, K, Maturity, 
%                 NoSteps, NoPaths)

%% INPUTS
%   S0     - Current price of the underlying asset.
%   r      - Annualized continuously compounded risk-free rate of return
%            over the life of the option, expressed as a positive decimal
%            number.
%   q      - Annualized continuously compounded yield rate
%
%   V0     - Current variance of the underlying asset
%   xi     - volatility of volatility
%   theta  - long-term mean
%   kappa  - rate of mean-reversion
%
%   strike      - Vector of strike prices of the option
%   Maturity    - Time to expiration of the option, expressed in years.
%   NoSteps     - Number of time steps per path
%   NoPaths     - Number of paths (Monte-Carlo simulations)


%% OUPUTS
%   call_prices     - Prices (i.e., value) of a vector of European call options.
%   std_err         - Standard deviation of the error due to the Monte-Carlo 
%                     simulation:
%                     (std_err = std(sample)/sqrt(length(sample)))

%% REFERENCES
% [AN06] Andersen, Leif. 2006. “Efficient Simulation of the Heston Stochastic Volatility Model 1 Introduction.”

%% ATTRIBUTION
% Christos  Delivorias
% The University of Edinburgh
% August 2012
