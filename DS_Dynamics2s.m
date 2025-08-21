%This code simulates the dynamical system on the generated network and
%outputs the dynamical variables as a function of sigma. The dynamical
%variables calculated here are V2 and W2 only, hence this code is specific
%to the study of the linear stability analysis. The main difference between
%the Dyn1 and Dyn2 codes is the initial condition Psi0 which is random in
%Dyn1 but in Dyn2 it's aligned with the selected eigenvector plus some
%noise. Also, in Dyn2s, Psi0 is initialised for every sigma.

tic

clear all

%%Initialize random seed
rng('shuffle');

%% STEP 1: DEFINE NETWORK
load('../Codes_Networks/RG_k12_N50')
addpath('../Codes_DS')

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

%%Initialise remaining parameters
Sigma = 0:0.01:15;
Tmax = 15; % max time to integrate
dt = 0.001; %time step
time = 0:dt:Tmax;
ratio = 2/3; %Order parameters will be annealed over T=[ratio*Tmax,Tmax]

W2 = zeros(length(Sigma),5001); % 5001 is the max value of idx3
V2 = zeros(length(Sigma),5001);

%% STEP 3: RUN DYNAMICS
for idx1 = 1:length(Sigma)
    
    sigma = Sigma(idx1);

    %%Initialise Psi for each new sigma
    rvect = randn(N+L,1); % a random vector
    rvect = rvect -  dot(rvect,seve(:,sel2))*seve(:,sel2); % remove random component along selected eigenstate as this contronlled for.
    nrvect = rvect/ norm(rvect,2); % normalise the random vector.
    epsilon = 0.1;
    Psi0 = (sqrt(1-epsilon^2)*seve(:,sel2)) + (epsilon*nrvect); % the coeffs of seve and nrvect have this form so that Psi0 is normalised.
    Psi = Psi0;

    Psi_dot_vect = zeros(N+L,length(time));
    Psi_vect = zeros(N+L,length(time));

    idx3 = 0;

    for idx2 = 1:length(time) %time loop

        %%RK4 method
        [k1{1}] = func(Psi, sigma, Deqn,Omega); %Eqs (32) in notes
        [k2{1}] = func(Psi+0.5*dt*k1{1},sigma,Deqn,Omega);
        [k3{1}] = func(Psi+0.5*dt*k2{1},sigma,Deqn,Omega);
        [k4{1}] = func(Psi+dt*k3{1},sigma,Deqn,Omega);
        Psi_dot = (k1{1}+2*k2{1}+2*k3{1}+k4{1})/6;
        del_Psi = dt*Psi_dot;
        %%Final update
        Psi = Psi + del_Psi;

        if idx2 >= length(time)*ratio
            idx3 = idx3+1;
            Psi_dot_vect(:,idx3) = Psi_dot;
            Psi_vect(:,idx3) = Psi;
        end

    end
    
    Psi_dot_vect = Psi_dot_vect(:,1:idx3);
    pPsi_dot_vect = seve'*Psi_dot_vect;
    pPsi_dot_vect(sel2,:) = []; % remove component along selected eigenstate
    V2(idx1,:) = mean((pPsi_dot_vect).^2);
    
    Psi_vect = Psi_vect(:,1:idx3);
    pPsi_vect = seve'*Psi_vect;
    pPsi_vect(sel2,:) = []; % remove component along selected eigenstate
    W2(idx1,:) = mean((pPsi_vect).^2);
end

%%Set name of output file
if net_type == "RG"
    formatSpec = sprintf('DS_Dyn2s_m%d_E%d_lam%d_%d_RG_k%d_N%d',round(m),round(E),ismember(E,slam),isequal(round(abs(E),3),abs(m)),c,N);
elseif net_type == "SBM"
    formatSpec = sprintf('DS_Dyn2s_m%d_E%d_lam%d_%d_SBM_M%d_N%d',round(m),round(E),ismember(E,slam),isequal(round(abs(E),3),abs(m)),M,N);
end
%%Save data to output file
save(formatSpec,"-v7.3");

toc
