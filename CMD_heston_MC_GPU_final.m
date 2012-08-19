function [Payoff, call,std_err, V, S ]=CMD_heston_GPU(S0, rho, V0, xi, theta, kappa, K, T, steps, paths, lambda, r, q,NAG)
% Time granulation
dt = T/steps ;

xmean(1:steps)=0;
sd(1:steps)=1;

Payoff=[];
call=[];
std_err=[];
V=[];
S=[];


% Assign each worker to a different GPU
spmd
gpuDevice(labindex);
end

% Pre-cache result arrays
S = parallel.gpu.GPUArray.zeros(paths,1) ;
S(:,1) = S0 ; % add the initial price
V = parallel.gpu.GPUArray.zeros(paths,1) ;
V(:,1) = V0 ; % add the initial var
C = parallel.gpu.GPUArray.zeros(paths,1);
Payoff = parallel.gpu.GPUArray.zeros(paths,1);



disp('Pre-cache arrays ready')
gamma1=1/2 ; % See [AN06] for more options
gamma2=1/2 ; % See [AN06] for more options 
 
% Pre-cache calculation constants
K0 = -rho*kappa*theta*dt/xi ;
K1 = gamma1*dt*(kappa*rho/xi-1/2)-rho/xi ;
K2 = gamma2*dt*(kappa*rho/xi-1/2)+rho/xi ;
K3 = gamma1*dt*(1-rho^2) ;
K4 = gamma2*dt*(1-rho^2) ;

disp('Kappas ready')
% Set the seed for the RNG in MATLAB


seed = 1234;
VA=[];
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
disp('RNG ready')

% Main Monte Carlo loop

for step = 1: steps
    Zn1 = parallel.gpu.GPUArray.randn(paths,1);
	Zn2 = parallel.gpu.GPUArray.randn(paths,1);
	Uv  = parallel.gpu.GPUArray.rand(paths,1);
    
    [S,V]=arrayfun(@myGPUFun, Zn1,Zn2,Uv,S,V,theta, kappa, dt, xi,r,K0,K1,K2,K3,K4);

end
disp('.')
Payoff = max(S-K,0);        % Call option
C = S - K;
call = mean(C);
Payoff = exp(-r*T)* mean(Payoff);
std_err = std(C)/sqrt(paths);
end

%% Internal function for the parallelisation on the GPU
%
function [S,V] = myGPUFun(Zn1,Zn2,Uv,S,V,theta, kappa, dt, xi,r,K0,K1,K2,K3,K4)
 Vt=V;
 
 A1=exp(-kappa*dt) ;
 
 m=theta+(Vt-theta)*A1 ;
 s_2=Vt*((xi^2*A1)/kappa)*(1-A1)+(theta*xi^2)*(1-A1)^2/(2*kappa) ;
 psi=s_2./(m.*m) ;
 psi_c=1.5 ; % As mentioned by Andersen
 IDX=(psi<=psi_c);
 
 b_2=2/psi-1+sqrt(2/psi)*sqrt(2/psi-1) ;
 a=m/(1+b_2) ;
 tmp=(sqrt(b_2)+Zn1);
 VTMP1=a*tmp*tmp ;
  
  p=(psi-1)./(psi+1) ;
  beta= (1-p)./m ;
  
  VTMP2=(log((1-p)./(1-Uv))./beta).*(Uv>p) ;
  V=VTMP1.*IDX+VTMP2.*(1-IDX);
  
 
 
 S = S .* exp(r*dt+ K0 + K1*Vt+K2*V + sqrt(K3*Vt + K4*V)*Zn2);
end


%% Internal Function for the variance
function [Vdt]=QEvariance(Vt,theta,kappa,dt,xi, Zn, Uv) 
% Descritise the variance using the QE Scheme by Andersen [AN06]
% Strongly reflect at 0. Ensures no negative variance values.
A1=exp(-kappa*dt) ;
A2=(xi^2*A1)/kappa ;
A3=1-A1 ;
m=theta+(Vt-theta)*A1 ;
s_2=Vt*A2*A3+(theta*xi^2)*A3^2/(2*kappa) ;
psi=s_2/(m*m) ;
psi_c=1.5 ; % As mentioned by Andersen

if psi<=psi_c
	b_2=2/psi-1+sqrt(2/psi)*sqrt(2/psi-1) ;
	a=m/(1+b_2) ;
    tmp=(sqrt(b_2)+Zn);
    Vdt=a*tmp*tmp ;
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
end