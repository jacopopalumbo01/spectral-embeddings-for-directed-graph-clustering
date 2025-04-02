clear; close all;
rng(42);

%% This script generates images and plots regarding Block-Acycles

% Grap properties
k         = 8;
n         = 100;
conn_prob = 0.4;

% Uniform block membership distribution
rho_uniform = ones(k,1);
rho_uniform = rho_uniform / sum(rho_uniform);

% Generate block connection probability
P = ConnectionProbabilityMatrix("acyclic", k, conn_prob);

% Generate the graph
[W, nodes] = GenerateGraph(n,k,rho_uniform,P,1);

% Plot the graph
PlotCyclic(W, k, transpose(nodes));

% Perform BAS
[cluster_index, ~] = BAS(W, k, "transition", 1, true, false);

% Generate the perturbed graph
[W, nodes] = GenerateGraph(n,k,rho_uniform,P,1,0.2);

% Plot the graph
PlotCyclic(W, k, transpose(nodes));

% Perform BAS
[cluster_index, ~] = BAS(W, k, "transition", 1, true, false);
