function [a,G,flag]=random_graph(N,c)
%Generate a Poisson network
% INPUTS 
% N number of nodes
% c average degree
% OUTPUTS
% a adjacency matrix
% G graph
% flag 1 if network is connected, 0 if network is not connected

x=rand(N);
p=c/(N-1); %probability of a link
x=x<=p;
a=triu(x,1);
a=sparse(a+a'); %adjacency

%Check the network is connected: both should be true.
G=graph(a);
[~,binsizes]=conncomp(G);
flag=0;
if (binsizes(1)==N)
    flag=1;
end

end
