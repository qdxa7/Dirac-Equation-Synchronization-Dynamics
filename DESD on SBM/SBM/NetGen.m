%This code generates the SBM network.


V = 50; %number of nodes per cluster
M = 4; %number of clusters
while((flag==0)&&(nrun<10))
    [A,G,flag]=block_model(V,M);
    nrun=nrun+1;
end

%%Set name of output file
formatSpec = sprintf('SBM_M%d_N%d',M,N);

%%Save data to output file
save(formatSpec,"-v7.3");