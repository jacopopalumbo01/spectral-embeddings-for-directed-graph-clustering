function clusters = SkewSymmetricClustering(W, k)
% SkewSymmetricClustering - Clustering method based on skew-symmetric
% decomposition of a digraph adjacency matrix.
%
% Inputs:
%   W : (n x n) Adjacency matrix of a digraph (can be weighted or unweighted)
%   k : Desired number of clusters
%
% Output:
%   clusters : (n x 1) Cluster assignment from k-means

    % Step 1: Construct the skew-symmetric matrix K = W - W'
    K = W - W';

    % Step 2: Determine l based on parity of k
    if mod(k, 2) == 0
        l = k;
    else
        l = k - 1;
    end

    % Step 3: Compute truncated SVD of K
    [U, ~, ~] = svds(K, l);

    % Step 4: k-means on rows of U
    clusters = kmeans(U, k, 'Replicates', 10, 'MaxIter', 300);
end
