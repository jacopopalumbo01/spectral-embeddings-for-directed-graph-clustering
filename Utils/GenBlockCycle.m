function [W, nodes] = GenBlockCycle(num_blocks, num_nodes, conn_prob)
    % Calculate probability of belonging to a block
    prob = 1.0 / num_blocks;
    
    % Divide the nodes into blocks
    nodes = floor(rand(num_nodes,1) ./prob) + 1;
    
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
            end
        end
    end
end