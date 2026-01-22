% This code calculates and plots the cumulative spectral density of the Hamiltonian operator (H)
% defined for a simple network of the user's choice with an NxN adjacency matrix A, where N is
% the number of nodes. The adjacency matrix A must be defined before running this code (e.g. by
% loading a .mat file containing A). This code can be used to visualise the smallest isolated eigenstates
% of H that satisfy E>|m| as those can induce topological cluster-synchronization in the DESD model.

% Inputs:
% A - NxN adjacency matrix of the network

% Outputs:
% A plot of the cumulative spectral density of H

% Parameters (can be modified by the user):
% m - mass constant

%%Define network
N = size(A,1); %number of nodes
L = nnz(A)/2; %number of links
G = graph(A); %generate graph object
B1 = incidence(G); %generate incidence matrix

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
energy = eig(H);

%%Sort the energies of H into non-negative (nnenergy) and negative (nenergy) energies:
[idx,~]=find(energy>=0);
nnenergy = energy(idx);
nnenergy=sort(nnenergy,'descend');
[~,~,nenergy]=find(energy.*(energy<0));
nenergy=sort(nenergy,'ascend');

%%Plot cumulative spectral density of H
figure();
hold on
plot(nnenergy,[1:numel(nnenergy)]/numel(energy),'o','Color','b','DisplayName',sprintf('$m=%0.1f$',m));
plot(nenergy,[1:numel(nenergy)]/numel(energy),'o','Color','b','HandleVisibility','off');
xlabel('$E$','Interpreter','Latex','FontSize',20,'FontWeight','bold')
ylabel('$\rho_c(E)$','Interpreter','Latex','FontSize',20,'FontWeight','bold')
set(gca,'FontSize',16,'TickLabelInterpreter','latex');
legend('Box','off','FontSize',16,'Interpreter','latex')
box on
hold off