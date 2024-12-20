clear all; clc;

addpaths_SC;
% Synthetically generate an unperturbed block-cycle
graph_name = "Unperturbed Block-Cycle graph";
fprintf("%s\n", graph_name);
num_blocks = 8; % number of blocks
num_nodes  = 100; % numbert of nodes
conn_prob  = 0.8; % connection probability between consecutive blocks

[W,nodes] = GenBlockCycle(num_blocks, num_nodes, conn_prob);

% Plot the graph
PlotCyclic(W, num_blocks, nodes, graph_name)

% Compute the transition probability matrix
P = TransitionMatrix(W);

% Plot sparsity pattern of W
figure;
spy(W);
title("Unperturbed Block Cycle Graph Adjacency Matrix");

% Get clusters
[cycle_eigvals, cycle_eigvecs] = BCS(W, P, num_blocks, true, true, graph_name);

% Extract the real and imaginary part 
% from the cycle eigenvectors
cycle_real = real(cycle_eigvecs);
cycle_imag = imag(cycle_eigvecs);
% The new data matrix is [num_nodes, 2xecycle_eigenvecs]
data_real_imag = [cycle_real, cycle_imag];

%% Step 3: K-means on the rows
% [tau, C] = LloydClusterSensitive(cycle_eigvecs, k, 500);
[cluster_index, centroids] = kmeans(data_real_imag, num_blocks, 'Distance', 'sqeuclidean');

%% Evaluation
fprintf("----------\nEvaluation\n----------\n")
[RCut, NCut, NMI, FScore] = ComputeMetrics(nodes,cluster_index,W);

fprintf("   RCut:       %f\n", RCut);
fprintf("   NCut:       %f\n", NCut);
fprintf("   NMI:        %f\n", NMI);
fprintf("   F-score:    %f\n", FScore);

fprintf("\n\n");

% Synthetically generate a perturbed block-cycle
graph_name = "Perturbed Block-Cycle graph";
fprintf("%s\n", graph_name);
num_blocks = 8; % number of blocks
num_nodes  = 100; % numbert of nodes
conn_prob  = 0.8; % connection probability between consecutive blocks
pert_prob  = 0.2; % perturbation probability


[W,nodes] = GenBlockCycle(num_blocks, num_nodes, conn_prob, true, pert_prob);

% Plot the graph
PlotCyclic(W, num_blocks, nodes, graph_name);

% Compute the transition probability matrix
P = TransitionMatrix(W);

% Plot sparsity pattern of W
figure;
spy(W);
title("Perturbed Block Cycle Graph Adjacency Matrix");

% Get clusters
[cycle_eigvals, cycle_eigvecs] = BCS(W, P, num_blocks, true, true, graph_name);

% Extract the real and imaginary part 
% from the cycle eigenvectors
cycle_real = real(cycle_eigvecs);
cycle_imag = imag(cycle_eigvecs);
% The new data matrix is [num_nodes, 2xecycle_eigenvecs]
data_real_imag = [cycle_real, cycle_imag];

%% Step 3: K-means on the rows
% [tau, C] = LloydClusterSensitive(cycle_eigvecs, k, 500);
[cluster_index, centroids] = kmeans(data_real_imag, num_blocks, 'Distance', 'sqeuclidean');

%% Evaluation
fprintf("----------\nEvaluation\n----------\n")
[RCut, NCut, NMI, FScore] = ComputeMetrics(nodes,cluster_index,W);

fprintf("   RCut:       %f\n", RCut);
fprintf("   NCut:       %f\n", NCut);
fprintf("   NMI:        %f\n", NMI);
fprintf("   F-score:    %f\n", FScore);
