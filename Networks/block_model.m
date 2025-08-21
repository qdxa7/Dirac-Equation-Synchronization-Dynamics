function [a,G,flag]=block_model(V,M)
%Generate a block model with M clusters
% INPUTS 
% V number of nodes per cluster
% M number of clusters
% OUTPUTS
% a adjacency matrix
% G graph
% flag 1 if network is connected, 0 if network is not connected

c1 = 15; %10; %8; %20; %average intra-cluster degree.
c2 = 1; %0.1; %average inter-cluster degree.
p1 = c1/(V-1); %probability of links within a cluster (this was originally p1=c1/V which I think isn't accurate).
p2 = c2/((M-1)*V); %probability of links between clusters (this was originally p2=c2/((M-1)*(V-1)) which I think isn't accurate).

for i = 0:(M-1)
    x1 = rand(V);
    a([1+i*V:(i+1)*V],[1+i*V:(i+1)*V]) = x1<p1;
    for j = i+1:M-1
        x1 = rand(V);
        a([1+i*V:(i+1)*V],[1+j*V:(j+1)*V]) = x1<p2;
        a([1+j*V:(j+1)*V],[1+i*V:(i+1)*V]) = x1<p2;
    end
end
a = triu(a,1);
a = sparse(a + a');

%Check the network is connected: both should be true.
G=graph(a);
[~,binsizes]=conncomp(G);
flag=0;
if (binsizes(1)==V*M)
    flag=1;
end


end

