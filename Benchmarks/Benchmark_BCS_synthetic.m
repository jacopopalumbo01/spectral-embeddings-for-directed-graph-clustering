% This script is used to benchmark the Block Cyclic Spectral clustering
% algorithm on synthetic block cyclic graphs

% Add paths and rng
clear all; close all; clc;
rng(42);
addpaths_SC;

% Write output to file
diary on;
diary('BCS_synthetic_output.txt');

fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');
fprintf('|                                                                     |\n');
fprintf('|  Benchmark of Block Cyclic Spectral clustering on synthetic graphs  |\n');
fprintf('|                                                                     |\n');
fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n\n');



%% Uniform Block Cyclic without perturbation
fprintf("Testing Uniform Block Cyclic without perturbation...\n")

% Path of matrices
folder_name = sprintf("%s/Data/Synthetic/UniformBlockCyclicNoPerturbation", pwd);

% Test dimensions
num_blocks = [3;5;8;10];
num_nodes  = [100;200;500;1000;2500;5000];
num_tests = size(num_nodes,1);

for i = 1:size(num_blocks,1)
    % Initialize metrics vectors
    NCut   = zeros(num_tests, 1);
    RCut   = zeros(num_tests, 1);
    NMI     = zeros(num_tests, 1);
    FScore = zeros(num_tests, 1);

    for j = 1:size(num_nodes,1)
        % Load adjacency matrix and labels
        W = load(sprintf("%s/%dblocks_%dnodes_adj.mat", folder_name, ...
            num_blocks(i,1), num_nodes(j,1))).W;
        nodes = load(sprintf("%s/%dblocks_%dnodes_lab.mat", folder_name, ...
            num_blocks(i,1), num_nodes(j,1))).labels;

        % Get clusters
        [cluster_index, ~] = BCS(W, num_blocks(i,1), ...
            false, false);
        
        % Compute and save metrics
        [RCut(j,1), NCut(j,1), NMI(j,1), FScore(j,1)] = ComputeMetrics( ...
            nodes,cluster_index,W);
    end

    % Print results
    fprintf("---------------------\n");
    fprintf("   %d blocks\n", num_blocks(i,1));
    fprintf("---------------------\n");
    % Generate table with computed metrics
    T = table(num_nodes, NCut, RCut, NMI, FScore); 
    % Display the table
    disp(T);
end


%% Uniform Block Cyclic with perturbation

fprintf("Testing Uniform Block Cyclic with perturbation...\n")

% Path of matrices
folder_name = sprintf("%s/Data/Synthetic/UniformBlockCyclicPerturbation", pwd);

% Test dimensions
num_blocks = [3;5;8;10];
num_nodes  = [2500];
pert_prob  = [0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.65; 0.7; 0.75; 0.8];
num_tests = size(pert_prob,1);

for i = 1:size(num_blocks,1)
    for j = 1:size(num_nodes,1)
        % Initialize metrics vectors
        NCut   = zeros(num_tests, 1);
        RCut   = zeros(num_tests, 1);
        NMI     = zeros(num_tests, 1);
        FScore = zeros(num_tests, 1);

        for z = 1:size(pert_prob,1)

        % Load adjacency matrix and labels
        W = load(sprintf("%s/%dblocks_%dnodes_%fpert_adj.mat", folder_name, ...
            num_blocks(i,1), num_nodes(j,1), pert_prob(z,1))).W;
        nodes = load(sprintf("%s/%dblocks_%dnodes_%fpert_lab.mat", folder_name, ...
            num_blocks(i,1), num_nodes(j,1), pert_prob(z,1))).labels;
    
        % Get clusters
        [cluster_index, ~] = BCS(W, num_blocks(i,1), ...
            false, false);
        
        % Compute and save metrics
        [RCut(z,1), NCut(z,1), NMI(z,1), FScore(z,1)] = ComputeMetrics( ...
            nodes,cluster_index,W);
        end
        % Print results
        fprintf("--------------------------\n");
        fprintf("   %d blocks, %d nodes    \n", num_blocks(i,1), num_nodes(j,1));
        fprintf("--------------------------\n");
        % Generate table with computed metrics
        T = table(pert_prob, NCut, RCut, NMI, FScore); 
        % Display the table
    disp(T);
    end
end

%% Block Cyclic without perturbation
fprintf("Testing Block Cyclic without perturbation...\n")

% Path of matrices
folder_name = sprintf("%s/Data/Synthetic/BlockCyclicNoPerturbation", pwd);

% Test dimensions
num_blocks = [3;5;8;10];
num_nodes  = [100;200;500;1000;2500;5000];
num_tests = size(num_nodes,1);

for i = 1:size(num_blocks,1)
    % Initialize metrics vectors
    NCut   = zeros(num_tests, 1);
    RCut   = zeros(num_tests, 1);
    NMI     = zeros(num_tests, 1);
    FScore = zeros(num_tests, 1);
    
    for j = 1:size(num_nodes,1)
        % Load adjacency matrix and labels
        W = load(sprintf("%s/%dblocks_%dnodes_adj.mat", folder_name, ...
            num_blocks(i,1), num_nodes(j,1))).W;
        nodes = load(sprintf("%s/%dblocks_%dnodes_lab.mat", folder_name, ...
            num_blocks(i,1), num_nodes(j,1))).labels;

        % Get cluster
        [cluster_index, ~] = BCS(W, num_blocks(i,1), ...
            false, false);
        
        % Compute and save metrics
        [RCut(j,1), NCut(j,1), NMI(j,1), FScore(j,1)] = ComputeMetrics( ...
            nodes,cluster_index,W);
    end

    % Print results
    fprintf("---------------------\n");
    fprintf("   %d blocks\n", num_blocks(i,1));
    fprintf("---------------------\n");

    % Print membership distribution
    distribution = load(sprintf("%s/%dblocks_distribution.mat", folder_name, ...
            num_blocks(i,1))).block_weights;
    fprintf("Membership distribution:")
    fprintf('[');
    fprintf('%g, ', distribution(1:end-1));
    fprintf('%g]\n', distribution(end));

    % Generate table with computed metrics
    T = table(num_nodes, NCut, RCut, NMI, FScore); 
    % Display the table
    disp(T);
end


%% Block Cyclic with perturbation

fprintf("Testing Block Cyclic with perturbation...\n")

% Path of matrices
folder_name = sprintf("%s/Data/Synthetic/BlockCyclicPerturbation", pwd);

% Test dimensions
num_blocks = [3;5;8;10];
num_nodes  = [2500];
pert_prob  = [0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.65; 0.7; 0.75; 0.8];
num_tests = size(pert_prob,1);

for i = 1:size(num_blocks,1)

    for j = 1:size(num_nodes,1)
        % Initialize metrics vectors
        NCut   = zeros(num_tests, 1);
        RCut   = zeros(num_tests, 1);
        NMI     = zeros(num_tests, 1);
        FScore = zeros(num_tests, 1);

        for z = 1:size(pert_prob,1)

        % Load adjacency matrix and labels
        W = load(sprintf("%s/%dblocks_%dnodes_%fpert_adj.mat", folder_name, ...
            num_blocks(i,1), num_nodes(j,1), pert_prob(z,1))).W;
        nodes = load(sprintf("%s/%dblocks_%dnodes_%fpert_lab.mat", folder_name, ...
            num_blocks(i,1), num_nodes(j,1), pert_prob(z,1))).labels;

        % Get clusters
        [cluster_index, ~] = BCS(W, num_blocks(i,1), ...
            false, false);
        
        % Compute and save metrics
        [RCut(z,1), NCut(z,1), NMI(z,1), FScore(z,1)] = ComputeMetrics( ...
            nodes,cluster_index,W);
        end
        % Print results
        fprintf("--------------------------\n");
        fprintf("   %d blocks, %d nodes    \n", num_blocks(i,1), num_nodes(j,1));
        fprintf("--------------------------\n");

        % Print membership distribution
        distribution = load(sprintf("%s/%dblocks_distribution.mat", folder_name, ...
            num_blocks(i,1))).block_weights;
        fprintf("Membership distribution:")
        fprintf('[');
        fprintf('%g, ', distribution(1:end-1));
        fprintf('%g]\n', distribution(end));

        % Generate table with computed metrics
        T = table(pert_prob, NCut, RCut, NMI, FScore); 
        % Display the table
    disp(T);
    end
end

% End output to file
diary off;