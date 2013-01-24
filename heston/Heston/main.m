%http://www.mathworks.com/matlabcentral/fileexchange/25771-heston-option-pricer


S0 = 100; r = 0.02;
    V0 = 0.04; eta = 0.7; theta = 0.06; kappa = 1.5;
    strike = 85:5:115; T = 0.25
    M = 2000; % Number of paths.
    N = 250; % Number of time steps per path

      [call_prices, std_errs] = Heston(S0, r, V0, eta, theta, kappa, strike, T, M, N)