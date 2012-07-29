%------------------------------------------------------------------%
function [ S ] = CIRpath( S0, mu, lambda, xi, T, dt, paths)
% CIRPATH    Generates paths based on the [ see \cite{Gillespie1996} 
%            for more on implementation]
%                                                                  %
%------------------------------------------------------------------%
% INPUTS     S0     - Initial value, positive scalar
%            mu     - Long term mean, positive scalar
%            lambda - Mean reversion rate, positive scalar
%            xi     - Volatility, positive scalar
%            T      - Maturity, positive scalar
%            dt     - Time intervals, positive scalar
%            paths  - Number of processes to create
%                                                                  %
%------------------------------------------------------------------%
% OUTPUTS    Plot of the CIR paths
%                                                                  %
%------------------------------------------------------------------%
% EXAMPLES   CIRpaths(5,1,3,0.6,10,1,5)
%                                                                  %
%------------------------------------------------------------------%
% Christos  Delivorias - OR MSc Student                            %
% The University of Edinburgh                                      %
% August 2012                                                      %
%------------------------------------------------------------------%
periods = floor(T /dt);

S = nan(paths,periods);

% Initial value
S(:,1) = S0;


lexpodt = exp(-lambda*dt);
for p=1:paths
    for t=2:periods
        S(p,t) = S(p,t-1)*lexpodt + mu*(1-lexpodt) + ...
            xi*sqrt((1-lexpodt^2)/2*lambda) * normrnd(0,1);
    end
end

plot(S','LineWidth',2);
