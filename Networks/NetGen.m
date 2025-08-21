%This code generates the network.

tic

clear all

%% STEP 1: DEFINE NETWORK
%%Generate network
net_type = "RG"; %choose network type

flag=0;
nrun=0;
if net_type == "RG" % fully connected
    N = 50; %number of nodes
    c = 12; %average degree
    while((flag==0)&&(nrun<10))
        [a,G,flag]=random_graph(N,c);
        nrun=nrun+1;
    end
elseif net_type == "SBM" %stochastic block model
    V = 50; %number of nodes per cluster
    M = 4; %number of clusters
    while((flag==0)&&(nrun<10))
        [a,G,flag]=block_model(V,M);
        nrun=nrun+1;
    end
    N = V*M; %total number of nodes
end

%spy(a);

%%Extract info about network
L = nnz(a)/2; %number of links
k = sum(a,2); %degree of each node
K = diag(k); %diagonal degree matrix
B1 = incidence(G); %incidence matrix
L0 = K - a; %graph Laplacian
L1d = B1'*B1; %first Hodge Laplacian
A1 = diag(diag(L1d)) - L1d; %adjacency matrix for links
kappa1 = sparse((abs(B1')*k)-2); %equal to k_{i}+k_{j}-2 for link [i,j]

%%Set name of output file
if net_type == "RG"
    formatSpec = sprintf('RG_k%d_N%d',c,N);
elseif net_type == "SBM"
    formatSpec = sprintf('SBM_M%d_N%d',M,N);
end

%%Save data to output file
save(formatSpec,"-v7.3");

toc