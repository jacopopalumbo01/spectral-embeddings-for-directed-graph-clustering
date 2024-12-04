clear all; clc;

addpaths_SC;
% Synthetically generate a block-cycle
num_blocks = 8; % number of blocks
num_nodes  = 100; % numbert of nodes
conn_prob  = 0.8; % connection probability between consecutive blocks

[W,nodes] = GenBlockCycle(num_blocks, num_nodes, conn_prob);

% Plot the graph
PlotCyclic(W, num_blocks, nodes);

% Compute the transition probability matrix
P = TransitionMatrix(W);

% Get clusters
[cycle_eigvals, cycle_eigvecs] = BCS(W, P, num_blocks);

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

%% Internal metrics (based only on graph structure and weights)
normalized = 1;
NCut = computeRCutValue(cluster_index,W,normalized);

RCut = computeRCutValue(cluster_index,W,~normalized);
fprintf("   RCut:       %f\n", RCut);
fprintf("   NCut:       %f\n", NCut);

%% External metrics (based on labels)
% Compute nmi
NMI = nmi(nodes, cluster_index);
fprintf("   NMI:        %f\n", NMI);

[inferred_labels,~] = label_data(cluster_index,nodes,1);

% Compute f-score
[Scores] = evaluate_scores(nodes,inferred_labels);
F_score = Scores(3);
fprintf("   F-score:    %f\n", F_score);