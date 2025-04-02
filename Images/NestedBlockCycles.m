clear; close all;
rng(42);

%% This script generates images and plots regarding Nested Block-Cycles

% Grap properties
k         = 4;
n         = 20;
conn_prob = 0.8;

% Uniform block membership distribution
rho_uniform = ones(k,1);
rho_uniform = rho_uniform / sum(rho_uniform);

% Generate block connection probability
P = ConnectionProbabilityMatrix("nested", k, conn_prob);

% Generate the graph
[W, nodes] = GenerateGraph(n,k,rho_uniform,P);

% Plot the graph
PlotCyclic(W, k, transpose(nodes));

% Perform BCS
[cluster_index, ~] = BAS(W, k, "transition", 1, true, false);

