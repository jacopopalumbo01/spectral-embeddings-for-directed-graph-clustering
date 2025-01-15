clear all; clc;

addpaths_SC;

% Synthetically generate an unperturbed uniform block-cycle
graph_name = "Unperturbed Uniform Block-Cycle graph";
fprintf("%s\n", graph_name);
num_blocks    = 8; % number of blocks
num_nodes     = 100; % numbert of nodes
conn_prob     = 0.8; % connection probability between consecutive blocks
block_weights = ones(num_blocks, 1) .* (1/num_blocks); % membership probability

[W,nodes] = GenBlockCycle(num_blocks, num_nodes, conn_prob, block_weights);

% Plot the graph
PlotCyclic(W, num_blocks, nodes, graph_name)

% Get clusters
[cluster_index, ~] = BCS(W, num_blocks, true, true, graph_name);


%% Evaluation
fprintf("----------\nEvaluation\n----------\n")
[RCut, NCut, NMI, FScore] = ComputeMetrics(nodes,cluster_index,W);

fprintf("   RCut:       %f\n", RCut);
fprintf("   NCut:       %f\n", NCut);
fprintf("   NMI:        %f\n", NMI);
fprintf("   F-score:    %f\n", FScore);

fprintf("\n\n");



% Synthetically generate an unperturbed uniform block-cycle
graph_name = "Unperturbed not Uniform Block-Cycle graph";
fprintf("%s\n", graph_name);
num_blocks    = 8; % number of blocks
num_nodes     = 100; % numbert of nodes
conn_prob     = 0.8; % connection probability between consecutive blocks
block_weights = transpose([0.18 0.2 0.05 0.2 0.14 0.04 0.07 0.12]); % membership probability

[W,nodes] = GenBlockCycle(num_blocks, num_nodes, conn_prob, block_weights);

% Plot the graph
PlotCyclic(W, num_blocks, nodes, graph_name)

% Get clusters
[cluster_index, ~] = BCS(W, num_blocks, true, true, graph_name);

%% Evaluation
fprintf("----------\nEvaluation\n----------\n")
[RCut, NCut, NMI, FScore] = ComputeMetrics(nodes,cluster_index,W);

fprintf("   RCut:       %f\n", RCut);
fprintf("   NCut:       %f\n", NCut);
fprintf("   NMI:        %f\n", NMI);
fprintf("   F-score:    %f\n", FScore);

fprintf("\n\n");


% Synthetically generate a perturbed block-cycle
graph_name = "Perturbed Uniform Block-Cycle graph";
fprintf("%s\n", graph_name);
num_blocks = 8; % number of blocks
num_nodes  = 100; % numbert of nodes
conn_prob  = 0.8; % connection probability between consecutive blocks
block_weights = ones(num_blocks, 1) .* (1/num_blocks); % membership probability
pert_prob  = 0.2; % perturbation probability


[W,nodes] = GenBlockCycle(num_blocks, num_nodes, conn_prob, block_weights, true, pert_prob);

% Plot the graph
PlotCyclic(W, num_blocks, nodes, graph_name);

% Get clusters
[cluster_index, ~] = BCS(W, num_blocks, true, true, graph_name);

%% Evaluation
fprintf("----------\nEvaluation\n----------\n")
[RCut, NCut, NMI, FScore] = ComputeMetrics(nodes,cluster_index,W);

fprintf("   RCut:       %f\n", RCut);
fprintf("   NCut:       %f\n", NCut);
fprintf("   NMI:        %f\n", NMI);
fprintf("   F-score:    %f\n", FScore);

fprintf("\n\n");


% Synthetically generate a perturbed block-cycle
graph_name = "Perturbed not Uniform Block-Cycle graph";
fprintf("%s\n", graph_name);
num_blocks = 8; % number of blocks
num_nodes  = 100; % numbert of nodes
conn_prob  = 0.8; % connection probability between consecutive blocks
block_weights = transpose([0.18 0.2 0.05 0.2 0.14 0.04 0.07 0.12]); % membership probability
pert_prob  = 0.2; % perturbation probability


[W,nodes] = GenBlockCycle(num_blocks, num_nodes, conn_prob, block_weights, true, pert_prob);

% Plot the graph
PlotCyclic(W, num_blocks, nodes, graph_name);

% Get clusters
[cluster_index, ~] = BCS(W, num_blocks, true, true, graph_name);

%% Evaluation
fprintf("----------\nEvaluation\n----------\n")
[RCut, NCut, NMI, FScore] = ComputeMetrics(nodes,cluster_index,W);

fprintf("   RCut:       %f\n", RCut);
fprintf("   NCut:       %f\n", NCut);
fprintf("   NMI:        %f\n", NMI);
fprintf("   F-score:    %f\n", FScore);

