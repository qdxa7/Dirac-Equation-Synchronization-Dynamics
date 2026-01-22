% This code simulates the DESD on a simple network of the user's choice with an NxN adjacency 
% matrix A, where N is the number of nodes. The code outputs the dynamical variables as a
% function of time. The adjacency matrix A must be defined before running this code (e.g. by
% loading a .mat file containing A).

% Inputs:
% A - NxN adjacency matrix of the network
% DESD.m - the auxiliary function that implements the vector field associated with DESD.

% Outputs:
% Psi_vs_time - phases of nodes and links as a function of time
% Psi_dot_vs_time - phase velocities of nodes and links as a function of time

% Parameters (can be modified by the user):
% m - mass constant
% sel - index of the eigenstate to be selected. Note that sel=1 selects the harmonic
%       eigenstate(s) with E=|m| and sel=1+n selects the nth isolated eigenstate with
%       E>|m| if the spectrum contains any. These isolated eigenstates can induce
%       topological cluster-synchronization in the DESD model. To know how many isolated
%       eigenstates exist in the spectrum of the Hamiltonian operator defined for your 
%       network of choice, run the Spectra.m code provided in this repository. Note also
%       that the non-negative energies are selected by default. To select the corresponding
%       negative ones, simply negate the variable E with a minus sign infront of 'energy(sel)'.
% Omega0, tau0 - mean and precision of the distribution of node intrinsic frequencies
% Omega1, tau1 - mean and precision of the distribution of link intrinsic frequencies
% cOmega(sel) - component of Omega in the direction of the selected eigenstate
% sigma - coupling constant in the DESD model
% Tmax - maximum time to integrate the dynamical system
% dt - time step for numerical integration

%% STEP 1: DEFINE NETWORK

N = size(A,1); %number of nodes
L = nnz(A)/2; %number of links
G = graph(A); %generate graph object
B1 = incidence(G); %generate incidence matrix

%% STEP 2: DEFINE DYNAMICAL SYSTEM

func = @DESD;

%%Generate the Topological Dirac operator (D):
D=sparse(N+L,N+L);
D(1:N, N+1:N+L)=B1;
D(N+1:N+L, 1:N)=transpose(B1);

%%Generate the gamma matrix (gamma):
gamma=sparse(N+L,N+L);
gamma(1:N, 1:N)=speye(N);
gamma(N+1:N+L, N+1:N+L)=-speye(L);

%%Set the mass constant (m):
m = 1;

%%Define the Hamiltonian operator (H):
H = D+(m*gamma);

%%Calculate the spectra of H:
[estate,energy] = eigs(H,N+L);
energy = real(diag(energy));
[energy,idx] = sort(energy,'ascend'); %sort energies from lowest to highest
estate = estate(:,idx); %sort eigenstates compatibly with energies
[idx,~]=find(energy>=0);
energy = energy(idx); %save only non-negative energies (negative energies are equal to non-negative energies with opposite sign)

%%Select desired eigenstate by setting 'sel' to the index of its associated energy (E): 
sel = 2;
E = energy(sel);

%%Define the final operator to be used in the DESD model (Deqn):
Deqn = H-(E*speye(N+L));

%%Define the intrinsic frequency vector (Omega):
Omega0 = 0; %this is the mean node intrinsic frequency
tau0 = 1; %this is the precision of the distribution of node intrinsic frequencies
w = Omega0 + randn(N,1)/tau0; %Gaussian node frequencies
Omega1 = 0; %this is the mean node intrinsic frequency
tau1 = 1; %this is the precision of the distribution of node intrinsic frequencies
wedge = Omega1 + randn(L,1)/tau1;
Omega = [w;wedge];

%%(Optional) set the component of Omega in the direction of the selected eigenstate to a specific value (e.g. 1):
cOmega = estate'*Omega; %The component of Omega in the direction of each eigenstate.
cOmega(sel) = 1;
Omega = (cOmega'*estate')';

%%Initialise random phase configuration (Psi0):
Psi0 = 2*pi*rand(N+L,1);
Psi = Psi0;

%%Set remaining DESD and simulation parameters:
sigma = 15;
Tmax = 30; % max time to integrate
dt = 0.001; %time step
time = 0:dt:Tmax;

%%Preallocate arrays to store dynamical variables:
Psi_dot_vs_time = zeros(N+L,length(time));
Psi_vs_time = zeros(N+L,length(time));

%% STEP 3: RUN DYNAMICS

for idx = 1:length(time) %time loop

    %%RK4 method
    [k1{1}] = func(Psi, sigma, Deqn,Omega);
    [k2{1}] = func(Psi+0.5*dt*k1{1},sigma,Deqn,Omega);
    [k3{1}] = func(Psi+0.5*dt*k2{1},sigma,Deqn,Omega);
    [k4{1}] = func(Psi+dt*k3{1},sigma,Deqn,Omega);
    
    Psi_dot = (k1{1}+2*k2{1}+2*k3{1}+k4{1})/6;
    Psi_dot_vs_time(:,idx) = Psi_dot;
    del_Psi = dt*Psi_dot;
    Psi = Psi + del_Psi;
    Psi_vs_time(:,idx) = Psi;

end

%%Set name of output file
fileName = sprintf('DESD_m%d_E%d',round(m),round(E));

%%Save data to output file
save(fileName,"-v7.3");