% This script is used to generate the adjacency matrix and the labels
% used in the evaluation of the spectral clustering methods.

% Add paths and rng
clear all; close all; clc;
rng(42);
addpaths_SC;

fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');
fprintf('|                                                           |\n');
fprintf('|  Generation of synthetic dataset for spectral clustering  |\n');
fprintf('|                                                           |\n');
fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n\n');

%%
fprintf("Generating synthetic uniform block-cyclic graph without perturbation\n");
folder_name = sprintf("%s/Data/Synthetic/UniformBlockCyclicNoPerturbation", pwd);

% Graph properties
num_blocks = [3;5;8;10];
num_nodes  = [100;200;500;1000;2500;5000];
conn_prob  = 0.8;

% Check and create folder where to store the data
if not(isfolder(folder_name))
   mkdir(folder_name)
end

% Generate data for each test case
for i = 1:size(num_blocks,1)
    for j = 1:size(num_nodes,1)
        % Compute uniform membership probability
        block_weights = ones(num_blocks(i,1), 1) .* (1/num_blocks(i,1));

        % Generate new graph
        [W, labels] = GenBlockCycle(num_blocks(i,1), num_nodes(j,1), ...
            conn_prob, block_weights);

        % Save graph adjacency and labels
        save(sprintf("%s/%dblocks_%dnodes_adj.mat", folder_name, ...
            num_blocks(i,1), num_nodes(j,1)), "W");

        save(sprintf("%s/%dblocks_%dnodes_lab.mat", folder_name, ...
            num_blocks(i,1), num_nodes(j,1)), "labels");

        fprintf("   Saved: %d blocks, %d nodes, %f conn_prob\n", ...
            num_blocks(i,1), num_nodes(j,1), conn_prob);
    end
end

fprintf("\nAll saved inside %s\n\n", folder_name);

%%
fprintf("Generating synthetic uniform block-cyclic graph with perturbation\n");
folder_name = sprintf("%s/Data/Synthetic/UniformBlockCyclicPerturbation", pwd);

% Graph properties
num_blocks = [3;5;8;10];
num_nodes  = [2500];
conn_prob  = 0.8;
pert_prob  = [0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.65; 0.7; 0.75; 0.8];

% Check and create folder where to store the data
if not(isfolder(folder_name))
   mkdir(folder_name)
end

% Generate data for each test case
for i = 1:size(num_blocks,1)
    for j = 1:size(num_nodes,1)
        for z = 1:size(pert_prob,1)
        % Compute uniform membership probability
        block_weights = ones(num_blocks(i,1), 1) .* (1/num_blocks(i,1));

        % Generate new graph
        [W, labels] = GenBlockCycle(num_blocks(i,1), num_nodes(j,1), ...
            conn_prob, block_weights, true, pert_prob(z,1));

        % Save graph adjacency and labels
        save(sprintf("%s/%dblocks_%dnodes_%fpert_adj.mat", folder_name, ...
            num_blocks(i,1), num_nodes(j,1), pert_prob(z,1)), "W");

        save(sprintf("%s/%dblocks_%dnodes_%fpert_lab.mat", folder_name, ...
            num_blocks(i,1), num_nodes(j,1), pert_prob(z,1)), "labels");

        fprintf("   Saved: %d blocks, %d nodes, %f conn_prob, %f perturbation\n", ...
            num_blocks(i,1), num_nodes(j,1), conn_prob, pert_prob(z,1));
    
        end
    end
end
fprintf("\nAll saved inside %s\n\n", folder_name);

%%
fprintf("Generating synthetic block-cyclic graph without perturbation\n");
folder_name = sprintf("%s/Data/Synthetic/BlockCyclicNoPerturbation", pwd);

% Graph properties
num_blocks = [3;5;8;10];
num_nodes  = [100;200;500;1000;2500;5000];

% Check and create folder where to store the data
if not(isfolder(folder_name))
   mkdir(folder_name)
end

% Generate data for each test case
for i = 1:size(num_blocks,1)
    % Compute random connection probability
    % First generate a random integer vector
    block_weights = randi([10,100],num_blocks(i,1),1);
    % Then normalize to make it a probability distribution
    block_weights = block_weights / sum(block_weights);
    % Save it
    save(sprintf("%s/%dblocks_distribution.mat", folder_name, num_blocks(i,1)) ...
        ,"block_weights");
    % Print the obtained distribution
    fprintf("%d blocks membership distribution:\n[", num_blocks(i,1));
    for t = 1:size(block_weights,1)
        fprintf("%f ", block_weights(t,1));
    end
    fprintf("]\n");

    for j = 1:size(num_nodes,1)
        % Generate new graph
        [W, labels] = GenBlockCycle(num_blocks(i,1), num_nodes(j,1), ...
            conn_prob, block_weights);

        % Save graph adjacency and labels
        save(sprintf("%s/%dblocks_%dnodes_adj.mat", folder_name, ...
            num_blocks(i,1), num_nodes(j,1)), "W");

        save(sprintf("%s/%dblocks_%dnodes_lab.mat", folder_name, ...
            num_blocks(i,1), num_nodes(j,1)), "labels");

        fprintf("   Saved: %d blocks, %d nodes, %f conn_prob\n", ...
            num_blocks(i,1), num_nodes(j,1), conn_prob);
    end
end

fprintf("\nAll saved inside %s\n\n", folder_name);


%%
fprintf("Generating synthetic block-cyclic graph with perturbation\n");
folder_name = sprintf("%s/Data/Synthetic/BlockCyclicPerturbation", pwd);

% Graph properties
num_blocks = [3;5;8;10];
num_nodes  = [2500];
pert_prob  = [0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.65; 0.7; 0.75; 0.8];

% Check and create folder where to store the data
if not(isfolder(folder_name))
   mkdir(folder_name)
end

% Generate data for each test case
for i = 1:size(num_blocks,1)
    % Compute random connection probability
    % First generate a random integer vector
    block_weights = randi([10,100],num_blocks(i,1),1);
    % Then normalize to make it a probability distribution
    block_weights = block_weights / sum(block_weights);
    % Save it
    save(sprintf("%s/%dblocks_distribution.mat", folder_name, num_blocks(i,1)) ...
        ,"block_weights");
    % Print the obtained distribution
    fprintf("%d blocks membership distribution:\n[", num_blocks(i,1));
    for t = 1:size(block_weights,1)
        fprintf("%f ", block_weights(t,1));
    end
    fprintf("]\n");

    for j = 1:size(num_nodes,1)
        for z = 1:size(pert_prob,1)
        % Generate new graph
        [W, labels] = GenBlockCycle(num_blocks(i,1), num_nodes(j,1), ...
            conn_prob, block_weights, true, pert_prob(z,1));

        % Save graph adjacency and labels
        save(sprintf("%s/%dblocks_%dnodes_%fpert_adj.mat", folder_name, ...
            num_blocks(i,1), num_nodes(j,1), pert_prob(z,1)), "W");

        save(sprintf("%s/%dblocks_%dnodes_%fpert_lab.mat", folder_name, ...
            num_blocks(i,1), num_nodes(j,1), pert_prob(z,1)), "labels");

        fprintf("   Saved: %d blocks, %d nodes, %f conn_prob, %f perturbation\n", ...
            num_blocks(i,1), num_nodes(j,1), conn_prob, pert_prob(z,1));
    
        end
    end
end

fprintf("\nAll saved inside %s\n\n", folder_name);
