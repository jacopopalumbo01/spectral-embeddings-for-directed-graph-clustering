clear all; clc;

% Synthetically generate a block-cycle
num_blocks = 8; % number of blocks
num_nodes  = 100; % numbert of nodes
conn_prob  = 0.8; % connection probability between consecutive blocks

[W,nodes] = GenBlockCycle(num_blocks, num_nodes, conn_prob);

% Plot the graph
PlotCyclic(W, num_blocks, nodes);

% Compute the transition matrix
P = TransitionMatrix(W);

% Get clusters
clusters = BCS(W, P, num_blocks);

%% Evaluation
fprintf("----------\nEvaluation\n----------\n")

%% Internal
% Compute RCut
RCut = computeRCutValue(clusters,W);
fprintf("   RCut:       %f\n", RCut);

%% External
% Compute nmi
NMI = nmi(nodes, clusters);
fprintf("   NMI:        %f\n", NMI);

[inferred_labels,~] = label_data(clusters,nodes,1);

% Compute f-score
[Scores] = evaluate_scores(nodes,inferred_labels);
F_score = Scores(3);
fprintf("   F-score:    %f\n", F_score);