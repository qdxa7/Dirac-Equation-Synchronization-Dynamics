%%The file SCmatrices88healthy.mat containing the data used was sourced from the paper "Human brain structural
% connectivity matrices–ready for modelling" (Nature, 2022). The data is publicly available at
% https://doi.org/10.17605/OSF.IO/YW5VF

%%Load structural connectivity (SC) matrix data
load('SCmatrices88healthy.mat')
%%The variable 'SCmatrices' is a 3D matrix where each slice SCmatrices(:,:,i)
%%corresponds to the structural connectivity matrix of a different healthy subject.
A_cell = squeeze(num2cell(SCmatrices, [2 3]));
%%Here, we select the SC matrix of the 10th subject as an example.
A = squeeze(A_cell{10});
A = (A+A')/2; % symmetrise
A = A>0.01; % select link probability threshold 
A = A - diag(diag(A)); % remove self-loops