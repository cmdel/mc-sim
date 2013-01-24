
function [call]=MCheston(So, r, V0, xi, theta, kappa, K, Maturity, NoSteps, NoPaths)

dt = Maturity/NoSteps;

% Initialise random variables


