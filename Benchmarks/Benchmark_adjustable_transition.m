clear;clc;
rng(47);

fprintf("------------------------------------------\n");
fprintf("BENCHMARK FOR ADJUSTABLE TRANSITION MATRIX\n");
fprintf("------------------------------------------\n");

%% Parameters
n         = 1000;         % Number of nodes
k         = 8;            % Number of blocks  
conn_prob = 0.7;

% Block membership distribution
rho = [0.18; 0.2; 0.05; 0.2; 0.14; 0.04; 0.07; 0.13];

% Betas to evaluate
betas = (0:0.1:2)';

%% Cyclic case
fprintf("Evaluating block-cyclic graph\n");
% Generate block connection probability
P = ConnectionProbabilityMatrix("cyclic", k, conn_prob);

% Generate graph
[W, nodes] = StochasticBlockmodel(n,k,rho,P);

num_experiments = size(betas,1);

modularities = zeros(num_experiments,1);

for i = 1:num_experiments
    [inferred_labels,~] = BCS(W, k, "power", betas(i));
    modularities(i) = Modularity(W,inferred_labels);
    
    fprintf("   beta=%f -> modularity=%f\n", betas(i), modularities(i));
end


fprintf("Evaluating block-acyclic graph\n");
% Generate block connection probability
P = ConnectionProbabilityMatrix("acyclic", k, conn_prob);

% Generate graph
[W, nodes] = StochasticBlockmodel(n,k,rho,P);

num_experiments = size(betas,1);

modularities = zeros(num_experiments,1);

for i = 1:num_experiments
    [inferred_labels,~] = BAS(W, k, "power", betas(i));
    modularities(i) = Modularity(W,inferred_labels);
    
    fprintf("   beta=%f -> modularity=%f\n", betas(i), modularities(i));
end