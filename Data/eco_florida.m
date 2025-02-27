clear;clc;
rng(47);

fprintf("------------------\n");
fprintf("EcoFlorida Dataset\n");
fprintf("------------------\n");

% Import the dataset
path = sprintf("%s/data/eco-florida/eco-florida.edges", pwd);
W = ImportEdges(path);

% Get the transition matrix
P = AcyclicTransitionMatrix(W);

% Plot the eigenvalues
[V, Cyclic_Eigs] = eigs(P,5,'lm');
Cyclic_Eigs = diag(Cyclic_Eigs);

D = eig(P);
modulus = abs(D);
    
% Select k largest eigenvalues and corresponding eigenvectors
% Note: We are interested in the eigenvalues with largest modulus
[~, indices] = maxk(modulus, 5);
cycle_eigvals = D(indices);
D(indices) = [];
cycle_eigvecs = V(:, indices);

PlotCyclicEig(D, cycle_eigvals, V, "");