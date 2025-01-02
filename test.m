clear all; clc;

addpaths_SC;

% This test is uniform because the membership probability is 1/num_blocks
% for each block (discrete uniform distribution).
fprintf("Testing synthetic uniform block-cyclic graph without perturbation\n");
num_blocks = [3;5;8;10];
num_nodes  = [100;200;500;1000;2500;5000];
conn_prob  = 0.8;

num_tests = size(num_nodes,1);

for i = 1:size(num_blocks,1)
    % Initialize metrics vectors
    NCut   = zeros(num_tests, 1);
    RCut   = zeros(num_tests, 1);
    NMI     = zeros(num_tests, 1);
    F_score = zeros(num_tests, 1);
    
    fprintf("")
    for j = 1:num_tests
        % Compute uniform membership probability
        block_weights = ones(num_blocks(i,1), 1) .* (1/num_blocks(i,1)); 
        % Generate new graph
        [W,nodes] = GenBlockCycle(num_blocks(i,1), num_nodes(j,1), conn_prob, block_weights);
    
        % Compute the transition matrix
        P = TransitionMatrix(W);
    
        % Get cycle eigenvalues
        [cycle_eigvals, cycle_eigvecs] = BCS(W, P, num_blocks(i,1), ...
            false, false);
        
        % Extract the real and imaginary part 
        % from the cycle eigenvectors
        cycle_real = real(cycle_eigvecs);
        cycle_imag = imag(cycle_eigvecs);
        % The new data matrix is [num_nodes, 2xecycle_eigenvecs]
        data_real_imag = [cycle_real, cycle_imag];

        % K-means on the rows
        [cluster_index, centroids] = kmeans(data_real_imag, num_blocks(i,1), 'Distance', 'sqeuclidean');

        % Compute and save metrics
        normalized = 1;
        NCut(j,1) = computeRCutValue(cluster_index,W,normalized);
        
        RCut(j,1) = computeRCutValue(cluster_index,W,~normalized);
        
        NMI(j, 1)   = nmi(nodes, cluster_index);
    
        [inferred_labels,~] = label_data(cluster_index,nodes,1);
        [Scores] = evaluate_scores(nodes,inferred_labels);
        F_score(j, 1) = Scores(3);
    end
    % Print results
    fprintf("---------------------\n");
    fprintf("   %d blocks\n", num_blocks(i,1));
    fprintf("---------------------\n");
    % Generate table with computed metrics
    T = table(num_nodes, NCut, RCut, NMI, F_score); 
    % Display the table
    disp(T);
end