%This code simulates the dynamical system on the generated network and
%outputs the dynamical variables as a function of time.

tic

clear all

%%Initialize random seed
rng('shuffle');

%% STEP 1: DEFINE NETWORK
load('../Codes_Networks/RG_k12_N1000.mat')


%% STEP 2: DEFINE DYNAMICAL SYSTEM
func = @Dirac;

%%Generate D and gamma and calculate spectra
D=sparse(N+L,N+L);
D(1:N, N+1:N+L)=B1;
D(N+1:N+L, 1:N)=transpose(B1);
m = 1;
gamma=sparse(N+L,N+L);
gamma(1:N, 1:N)=speye(N);
gamma(N+1:N+L, N+1:N+L)=-speye(L);
[eve,eva] = eigs((D+(m*gamma)),N+L);
lambda = real(diag(eva));
[slam,idx] = sort(lambda,'ascend'); %sort eigenvalues from lowest to highest
seve = eve(:,idx); %sort eigenvectors compatibly with eigenvalues
[~,~,neglam]=find(slam.*(slam<0));
[idx,~]=find(slam>=0);
nneglam = slam(idx);

% [~,sel1] = min(abs(nneglam - (pi*8/7))); %sel for selected
sel1 = 2;
E = nneglam(sel1);
sel2 = find(slam==E);
Deqn = D+(m*gamma)-(E*speye(N+L));

%%Define Omega
Omega0 = 0;
tau0 = 1;
w = Omega0 + randn(N,1)/tau0; %Gaussian node frequencies
Omega1 = 0;
tau1 = 1;
wedge = Omega1 + randn(L,1)/tau1;
Omega = ([w;wedge]'*seve')';

%%Design Omega
pOmega = seve'*Omega; %The component of Omega in the direction of each e.vect.
pOmega(sel2) = 1;
Omega = (pOmega'*seve')';

%%Initialise Psi
Psi0 = 2*pi*rand(N+L,1); %eve(:,sel); %2*pi*ones(N+L,1);
Psi = Psi0;

%%Initialise remaining parameters
sigma = 15;
Tmax = 30; % max time to integrate
dt = 0.001; %time step
time = 0:dt:Tmax;
ratio = 1/length(time); %Order parameters will be annealed over T=[ratio*Tmax,Tmax]

Psi_dot_tot = zeros(N+L,length(time));
del_Psi_tot = zeros(N+L,length(time));
Psi_tot = zeros(N+L,length(time));
theta = zeros(N,length(time));
phi = zeros(L,length(time));

idx1 = 0;


%% STEP 3: RUN DYNAMICS
for idx2 = 1:length(time) %time loop

    %%RK4 method
    [k1{1}] = func(Psi, sigma, Deqn,Omega); %Eqs (32) in notes
    [k2{1}] = func(Psi+0.5*dt*k1{1},sigma,Deqn,Omega);
    [k3{1}] = func(Psi+0.5*dt*k2{1},sigma,Deqn,Omega);
    [k4{1}] = func(Psi+dt*k3{1},sigma,Deqn,Omega);

    Psi_dot = (k1{1}+2*k2{1}+2*k3{1}+k4{1})/6;
    Psi_dot_tot(:,idx2) = Psi_dot;
    del_Psi = dt*Psi_dot;
    del_Psi_tot(:,idx2) = del_Psi;

    %%Final update
    Psi = Psi + del_Psi;
    Psi_tot(:,idx2) = Psi;
    theta(:,idx2) = Psi(1:N);
    phi(:,idx2) = Psi(N+1:N+L);

end

%%Set name of output file
if net_type == "RG"
    formatSpec = sprintf('DS_Dyn1t_m%d_E%d_lam%d_%d_RG_k%d_N%d',round(m),round(E),ismember(E,slam),isequal(round(abs(E),3),abs(m)),c,N);
elseif net_type == "SBM"
    formatSpec = sprintf('DS_Dyn1t_m%d_E%d_lam%d_%d_SBM_M%d_N%d',round(m),round(E),ismember(E,slam),isequal(round(abs(E),3),abs(m)),M,N);
end
%%Save data to output file
save(formatSpec,"-v7.3");

toc