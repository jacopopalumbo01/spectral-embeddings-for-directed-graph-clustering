function [W, nodes] = GenBlockCycle(num_blocks, num_nodes, conn_prob, block_weights, perturbed, pert_prob)
% GenBlockCycle - Generates a Block Cyclic graph
%% Syntax:
%        GenBlockCycle(num_blocks, num_nodes, conn_prob)
%
%% Input Arguments:
%       *Required Input Arguments*
%       - num_blocks:       Number of blocks of the Block Cyclic graph 
%       - num_nodes:        Number of nodes of the Block Cyclic graph
%       - conn_prob:        Connection probability between two consecutive
%                           blocks
%       - block_weights:    Vector containing weight for each block
%                           (probability for a node to belonging to each
%                           block). sum(block_weights) should be equal to 1
%       - perturbed:        Generate a perturbed Block Cycle graph
%                           (true/false)
%       - pert_prob:        Probability used in the perturbation phase to
%                           determine if an edge should be added or not. A
%                           value in [0,1]. 0 corresponds to no
%                           perturbation and 1 corresponds to an adjacency
%                           matrix with all cells set to 1.
%                           
%% Output:
%       - W:                Adjacency matrix
%       - nodes             Nodes membership
%

% Default case: no perturbation
if nargin < 5
    perturbed = false;
end

% Default perturbed case: 50% perturbation probability
if nargin < 6
    pert_prob = 0.5;
end

% Check if block_weights is valid
if round(sum(block_weights), 2) ~= 1.0
    error("Expected sum(block_weights) = 1.0, got %f", sum(block_weights));
else if size(block_weights, 1) ~= num_blocks
    error("Mismatch in dimensions: num_blocks = %d, " + ...
        "while size(block_weights, 1) = %d", num_blocks, size(block_weights, 1));
end


% Divide the nodes into blocks
nodes = randsample(num_blocks, num_nodes, true, block_weights);

% Order the nodes (i.e. all the first block at the start, followed by
% the second block and so on)
%nodes = sort(nodes);

% Initialize the adjacency matrix
W = zeros(num_nodes, num_nodes);

% Populate the adjacency matrix
for i = 1:num_nodes
    % Get block of the node
    block = nodes(i);
    
    % Calculate the consecutive block
    if block == num_blocks % The last block has the first block as consecutive
        next_block = 1;
    else
        next_block = block + 1;
    end

    % Connect with next block
    for j = 1:num_nodes
        if nodes(j) == next_block && rand(1) <= conn_prob
            W(i,j) = 1;
        else

            % Add perturbation
            if perturbed
                if rand(1) <= pert_prob
                    W(i,j) = 1;
                end
            end
        end
    end
end


end