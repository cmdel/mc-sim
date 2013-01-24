%% Monte Carlo Simulation of the Heston Model
% Calculate the value of a European call vanilla option
% using the Heston model for stochastic volatility

%% USAGE
% [payoff, call,std_err,V,S]=CMD_heston_MC(S0, rho, V0, xi, theta, kappa, K, T, 
% steps, paths, lambda, r, q, NAG)

%% INPUTS
%   S0     - Current price of the underlying asset.
%   r      - Annualized continuously compounded risk-free rate of return over the life of the option, expressed as a positive decimal number.
%   q      - Annualized continuously compounded yield rate
%
%   V0     - Current variance of the underlying asset
%   xi     - volatility of volatility
%   theta  - long-term mean
%   kappa  - rate of mean-reversion
%
%   K	   - Strike price of the option at maturity
%   T      - Time to expiration of the option, expressed in years.
%   steps  - Number of time steps per path
%   paths  - Number of paths (Monte-Carlo simulations)
%   lambda - The market risk aversion cost


%% OUPUTS
%   Payoffs      - Payoff prices given the strike price. max(S-K,0)
%   call_prices  - Prices (i.e., value) of a vector of European call options.
%   std_err      - Standard deviation of the error due to the Monte-Carlo simulation: (std_err = std(sample)/sqrt(length(sample)))
%   V            - The variance of the price of the underlying of all the paths
%   S            - The prices of the underlying for all the paths

%% REFERENCES
% [AN06] Andersen, Leif. 2006. “Efficient Simulation of the Heston 
% Stochastic Volatility Model.”

%% ATTRIBUTION
% Christos  Delivorias
% The University of Edinburgh
% August 2012
