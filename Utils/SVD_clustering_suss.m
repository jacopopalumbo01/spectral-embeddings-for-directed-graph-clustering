function [cluster_indexs, centroids] = SVD_clustering_suss(A, k)
% Code for Clustering Using SVD as per: https://arxiv.org/pdf/1108.2228
% Input: A adjacency matrix
%        k num_blocks

% Step 1: Compute the Singular Value Decomposition
[U, S, V] = svd(A);

% Step 2: Reduce Dimensionality
% Select the top-k singular vectors
U_k = U(:, 1:k);   % Left singular vectors
V_k = V(:, 1:k);   % Right singular vectors
S_k = S(1:k, 1:k); % Top-k singular values

% Coordinate scaled singular vector matrix
Z = [U_k * sqrt(S_k), V_k * sqrt(S_k)];

% Step 3: Cluster Z using kmeans

% Apply k-means clustering to the embeddings
[cluster_indexs, centroids] = kmeans(Z, k,'Replicates', 20);

end