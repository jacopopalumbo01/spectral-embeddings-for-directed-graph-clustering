function [k] = EstimateNumBlocksAcyclicWithModularity(W, max_nclust)
% EstimateNumBlocksAcyclicWithModularity - Compute number of clusters given 
% the adjacency matrix using the modularity as metric.
%
% Input:
%   - A:            Adjacency matrix
%   - max_k         Maximum number of clusters to consider
%
% Output:
%   - k:            Estimated number of clusters

modularities = zeros(max_nclust,1);

for k=2:max_nclust
    %% Step 1: Perform BAS for each number of clusters
    [cluster_index, ~] = BAS(W, k, false, false);
    
    %% Step 2: Compute modularity
    modularities(k - 1) = Modularity(W,cluster_index);
end

%% Step 4: Select candidate k with maximum modulus modularity 
[~, optimal_k_idx] = max(abs(modularities));
k = optimal_k_idx + 1;  % because candidates start at 2

end
