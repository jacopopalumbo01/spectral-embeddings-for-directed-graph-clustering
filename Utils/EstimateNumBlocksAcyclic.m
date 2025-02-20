function [k] = EstimateNumBlocksCyclic(A, max_k)
% EstimateNumBlocks - Compute number of clusters given the adjacency matrix.
%
% Input:
%   - A:            Adjacency matrix
%   - max_k         Maximum number of clusters to consider
%
% Output:
%   - k:            Estimated number of clusters

%% Step 1: Compute the transition matrix
P = AcyclicTransitionMatrix(A);

%% Step 2: Compute eigenvalues and eigenvectors
[eigvecs, D] = eig(P); 
eigvals = diag(D);

% Apply the filtering criteria
valid_indices = (real(eigvals) < 1) & (imag(eigvals) >= 0);
eigvals = eigvals(valid_indices);
eigvecs = eigvecs(:, valid_indices);

% Sort eigenvalues by descending absolute value so that the largest (in modulus) come first
[~, sort_idx] = sort(abs(eigvals), 'descend');
eigvals = eigvals(sort_idx);
eigvecs = eigvecs(:, sort_idx);

% Keep only the first k_to_find eigenvalues/eigenvectors.
% Here we define k_to_find as the smaller of floor(max_k/2) and the available number.
k_to_find = min(floor(max_k/2), length(eigvals));
% (Ensure at least 2 for a valid embedding.)
k_to_find = max(2, k_to_find);

% Truncate the eigenpairs
eigvals = eigvals(1:k_to_find);
eigvecs = eigvecs(:, 1:k_to_find);

%% Step 3: Evaluate candidate clusterings using the silhouette score
% We test candidate cluster numbers from 2 up to max_k
% (Note: if max_k exceeds 2*k_to_find, the embedding dimension remains constant.)
silhouette_scores = zeros(max_k - 1, 1);

for i = 2:max_k
    % Use the first floor(i/2) eigenvectors if available.
    num_cols = min(floor(i/2), k_to_find);
    submat = eigvecs(:, 1:num_cols);
    
    % Create a real embedding by concatenating real and imaginary parts
    submat_real_imag = [real(submat), imag(submat)];
    
    % Cluster into i clusters using k-means (with multiple replicates for stability)
    [clusters, ~] = kmeans(submat_real_imag, i, 'Distance', 'sqeuclidean', 'Replicates', 20, 'MaxIter', 10000);
    
    % Preallocate arrays for silhouette components
    eta = zeros(size(submat,1),1);   % average distance to other points in same cluster
    theta = zeros(size(submat,1),1); % minimum distance to points in other clusters
    s_values = zeros(size(submat,1),1); % silhouette value for each point
    
    for j = 1:size(submat,1)
        % Compute distances from point j to all other points in the embedding space
        point = [real(submat(j,:)), imag(submat(j,:))];
        embedding = [real(submat), imag(submat)];
        dists = pdist2(point, embedding);
        
        % Remove self-distance
        dists(j) = [];
        temp_clusters = clusters;
        temp_clusters(j) = [];
        
        % Within-cluster distances for point j
        same_cluster = dists(temp_clusters == clusters(j));
        if isempty(same_cluster)
            eta(j) = 0;  % singleton case
        else
            eta(j) = mean(same_cluster);
        end
        
        % Distances to points in other clusters
        other_clusters = dists(temp_clusters ~= clusters(j));
        if isempty(other_clusters)
            theta(j) = eta(j); % avoid NaN
        else
            theta(j) = min(other_clusters);
        end
        
        % Compute silhouette for point j (handle division by zero)
        if max(eta(j), theta(j)) == 0
            s_values(j) = 0;
        else
            s_values(j) = (theta(j) - eta(j)) / max(theta(j), eta(j));
        end
    end
    
    % Average silhouette score for candidate i
    silhouette_scores(i-1) = mean(s_values);
end

%% Step 4: Select candidate k with maximum silhouette score
[~, optimal_k_idx] = max(silhouette_scores);
k = optimal_k_idx + 1;  % because candidates start at 2

end
