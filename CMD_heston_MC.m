function [call,std_err, V, S, Payoff]=CMD_heston_MC(S0, rho, V0, xi, theta, kappa, K, T, steps, paths, lambda, r, q, NAG,Quasi)
% Time granulation
dt = T/steps ;
genid = 3;  % 1,2 = Sobol, 3 = Neidereitter, 4 = Faure
%Z = zeros(2*steps,1);
%Uv = zeros(2*steps,1);
% Memory pre-caching
% Pre-cache Normal(0,1) variates
if NAG
		if Quasi
			[Z, iref, ifail] = CMD_Normal_Quasi_NAG(0,1,2*steps,paths,genid) ;  % Create the N(0,1) with NAG
			Z = Z';
			Uv = CMD_Uniform_Quasi_NAG(steps, genid) ;
		else
			% TODO fix the normal NAG variates
			for p = 1:paths
				Z = [Z; CMD_Normal_NAG(0,1,2*steps,genid)];
				Uv =[Uv; CMD_Uniform_NAG(steps,genid)]; 	
			end
		end
	else
		Z = normrnd(0,1,paths,2*steps)	;				 % Create N(0,1) with Matlab
		Uv = rand(1,steps);
	end
size(Z);
Zn1 = Z(:,1:steps);			% Split the variates for the two Zn that are required
Zn2 = Z(:,steps+1:end);     % Split the variates for the two Zn that are required
size(Zn1);
size(Zn2);
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

% Main Monte Carlo loop
for pth = 1: paths
    for ts = 1:steps
		V(pth,ts+1) = QEvariance(V(pth,ts), theta, kappa, dt, xi, Zn1(pth,ts), Uv(ts)) ;
		St = S(pth,ts);
		Vt = V(pth,ts);
		Vdt = V(pth,ts+1);
		S(pth,ts+1) = St * exp( K0 + K1*Vt) * exp(K2*Vdt + sqrt(K3*Vt + K4*Vdt)*Zn2(pth,ts)); 
    end
	Payoff(pth) = max(S(pth,end)-K,0);
	C(pth) = S(pth,end);
end
call = mean(C);
Payoff = mean(Payoff);
std_err = std(C)/sqrt(paths);



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
