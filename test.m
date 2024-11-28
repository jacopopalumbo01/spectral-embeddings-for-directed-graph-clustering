clear all; clc;

num_blocks = 4;
num_nodes  = [10;20;50;100;200;500;1000];
conn_prob  = 0.8;

num_tests = size(num_nodes,1);

% Initialize metrics vectors
R_cut   = zeros(num_tests, 1);
NMI     = zeros(num_tests, 1);
F_score = zeros(num_tests, 1);

for i = 1:num_tests
    % Generate new graph
    [W,nodes] = GenBlockCycle(num_blocks, num_nodes(i,1), conn_prob);

    % Compute the transition matrix
    P = TransitionMatrix(W);

    % Get clusters
    clusters = BCS(W, P, num_blocks, false);
    
    % Compute and save metrics
    R_cut(i, 1) = computeRCutValue(clusters,W);
    NMI(i, 1)   = nmi(nodes, clusters);

    [inferred_labels,~] = label_data(clusters,nodes,1);
    [Scores] = evaluate_scores(nodes,inferred_labels);
    
    F_score(i, 1) = Scores(3);
end

figure;
semilogx(num_nodes,R_cut, "r-x");
xlabel("Number of nodes");
ylabel("Ratio Cut");
title("Ratio Cut")

figure;
semilogx(num_nodes,NMI, "r-x");
xlabel("Number of nodes");
ylabel("NMI");
title("NMI");

figure;
semilogx(num_nodes,NMI, "r-x");
xlabel("Number of nodes");
ylabel("F-Score");
title("F-score")