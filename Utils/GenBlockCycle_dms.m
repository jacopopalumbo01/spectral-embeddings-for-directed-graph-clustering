function [W, labels] = GenBlockCycle_dms(num_blocks, num_nodes, conn_prob, block_weights, perturbed, pert_prob)
    % GenBlockCycle - Generates a Block Cyclic graph
    %
    % Input:
    %   - num_blocks:  Number of blocks
    %   - num_nodes:   Total number of nodes
    %   - conn_prob:   Probability of connections between consecutive blocks
    %   - block_weights: Probability distribution for block membership
    %   - perturbed:   Boolean to add random perturbations
    %   - pert_prob:   Probability of random edge perturbations
    %
    % Output:
    %   - W: Adjacency matrix (NxN)
    %   - labels: Block membership labels for each node
    
        % Defaults
        if nargin < 5, perturbed = false; end
        if nargin < 6, pert_prob = 0.0; end
    
        % Validate block_weights
        if round(sum(block_weights), 2) ~= 1.0
            error("Sum of block_weights must be 1.0, got %f instead", sum(block_weights));
        elseif length(block_weights) ~= num_blocks
            error("block_weights length must match num_blocks.");
        end
    
        % Assign nodes to blocks based on block_weights
        labels = randsample(num_blocks, num_nodes, true, block_weights);
    
        % Initialize adjacency matrix
        W = zeros(num_nodes, num_nodes);
    
        % Populate adjacency matrix with cyclic block connections
        for i = 1:num_nodes
            current_block = labels(i);                       % Current node's block
            next_block = mod(current_block, num_blocks) + 1; % Cyclically get the next block
    
            for j = 1:num_nodes
                if labels(j) == next_block && rand() <= conn_prob
                    W(i, j) = 1; % Connect to nodes in the next block
                end
            end
        end
    
        % Add perturbations if enabled
        if perturbed
            random_edges = rand(num_nodes, num_nodes) <= pert_prob; % Random edge probability
            W = W | random_edges;                                  % Add perturbations
            W(1:num_nodes+1:end) = 0;                              % Remove self-loops (if any)
        end
    end
    