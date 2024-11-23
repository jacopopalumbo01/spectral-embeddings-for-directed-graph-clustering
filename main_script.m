clear all; clc;

% Synthetically generate a block-cycle
num_blocks = 3; % number of blocks
num_nodes  = 9; % numbert of nodes
conn_prob  = 0.8; % connection probability between consecutive blocks

[W,nodes] = GenBlockCycle(num_blocks, num_nodes, conn_prob);

% Plot the graph
PlotCyclic(W, num_blocks, nodes);

% Compute the transition matrix
P = TransitionMatrix(W);

% Compute the eigenvalues of the transition matrix
[V, D] = eig(P);

% Plot the cycle eigenvalues
figure;
scatter(real(D), imag(D), "rx");

% Get clusters
clusters = BCS(W, P, num_blocks);