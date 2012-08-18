function [Payoff, call,std_err, V, S ]=CMD_heston_MC(S0, rho, V0, xi, theta, kappa, K, T, steps, paths, lambda, r, q,NAG)
% Time granulation
dt = T/steps ;

xmean(1:steps)=0;
sd(1:steps)=1;

% Status output
disp('.')

% Pre-cache result arrays
S = zeros(paths,steps+1) ;
S(:,1) = S0 ; % add the initial price
V = zeros(paths,steps+1) ;
V(:,1) = V0 ; % add the initial var
C = zeros(paths,1);
Payoff = zeros(paths,1);

gamma1=1/2 ; % See [AN06] for more options
gamma2=1/2 ; % See [AN06] for more options 
 
% Pre-cache calculation constants
K0 = -rho*kappa*theta*dt/xi ;
K1 = gamma1*dt*(kappa*rho/xi-1/2)-rho/xi ;
K2 = gamma2*dt*(kappa*rho/xi-1/2)+rho/xi ;
K3 = gamma1*dt*(1-rho^2) ;
K4 = gamma2*dt*(1-rho^2) ;


if (NAG==1 || NAG==2)
    seed = [int64(1762543)];
    % genid and subid identify the base generator
    genid = int64(1);
    % genid = 3;  % 1,2 = Sobol, 3 = Neidereitter, 4 = Faure
    % Init the Quasi-random NAG generator for Normal(0,1) RVs.
    subid =  int64(1);
    stype = int64(1);   % Owen type scrambling of QRVs. 
						% 2 for Faure-Tezuka(FT), 3 for Owen and FT  
    idim = steps*3;
    iskip = int64(1000);
    nsdigi = int64(0);
    [state, ifail] = g05kf(genid, subid, seed);
    if(NAG==1)
      [iref, state, ifail] = g05yn(genid, stype, int64(idim), iskip, nsdigi, state);
    else
      [iref,ifail]=g05yl(int64(genid), int64(idim),int64(1));
    end
	[VA, iref,ifail] = g05ym(int64(paths),int64(idim),iref); % Create 3 URV
    VA=VA';
end

% Main Monte Carlo loop
for pth = 1: paths
	if (NAG==1 || NAG==2)
         Uv = VA(1:steps,pth);
         Zn1 = sqrt(2)*erfinv(2*VA(steps+1:2*steps,pth)-1); % Calculate N(0,1) with inverse CDF
         Zn2 = sqrt(2)*erfinv(2*VA(steps*2+1:end,pth)-1); % Calculate N(0,1) with inverse CDF
    elseif(NAG==-1)
        Zn1=randn(1,steps);
        Zn2=randn(2,steps);
		Uv=rand(1,steps);
	else
		Zn1=randn(1,steps/2);
		Zn1=[Zn1 -Zn1];
		Zn2=randn(1,steps/2);
		Zn2=[Zn2 -Zn2];
		Uv=rand(1,steps);
	end
    for ts = 1:steps
		V(pth,ts+1) = QEvariance(V(pth,ts), theta, kappa, dt, xi, Zn1(ts), Uv(ts));
		St = S(pth,ts);
		Vt = V(pth,ts);
		Vdt = V(pth,ts+1);
		S(pth,ts+1) = St * exp(r*dt+ K0 + K1*Vt) * exp(K2*Vdt + sqrt(K3*Vt + K4*Vdt)*Zn2(ts)); 
    end
	Payoff(pth) = max(S(pth,end)-K,0);        % Call option
	C(pth) = S(pth,end) - K;
end

call = mean(C);
Payoff = exp(-r*T)* mean(Payoff);
%std_err = std(C)/sqrt(paths);
for ts = 1:steps
	std_err(ts)=std(S(:,ts))/sqrt(paths);
end




% Internal Function for the variance
function [Vdt]=QEvariance(Vt,theta,kappa,dt,xi, Zn, Uv) 
% Descritise the variance using the QE Scheme by Andersen [AN06]
% Strongly reflect at 0. Ensures no negative variance values.
A1=exp(-kappa*dt) ;
A2=(xi^2*A1)/kappa ;
A3=1-A1 ;
m=theta+(Vt-theta)*A1 ;
s_2=Vt*A2*A3+(theta*xi^2)*A3^2/(2*kappa) ;
psi=s_2/m^2 ;
psi_c=1.5 ; % As mentioned by Andersen

if psi<=psi_c
	b_2=2/psi-1+sqrt(2/psi)*sqrt(2/psi-1) ;
	a=m/(1+b_2) ;
    Vdt=a*(sqrt(b_2)+Zn)^2 ;
else
	p=(psi-1)/(psi+1) ;
	beta= (1-p)/m ;
	if Uv<=p
		PSIu=0;
	else
		PSIu=log((1-p)/(1-Uv))/beta ;
	end
	Vdt=PSIu ;
end
