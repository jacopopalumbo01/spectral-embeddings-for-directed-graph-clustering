function [cluster_indexs, centroids] = SVD_clustering(A, k)
% Code for Clustering Using SVD
% Input: A adjacency matrix
%        k num_blocks

% Step 1: Compute the Singular Value Decomposition
[U, S, V] = svd(A);

% Step 2: Reduce Dimensionality
% Select the top-k singular vectors
U_k = U(:, 1:k); % Left singular vectors
V_k = V(:, 1:k); % Right singular vectors

% Step 3: Perform Clustering
% Combine embeddings for clustering (either U_k or V_k can be used depending on the task)
embeddings = U_k; % Use left singular vectors as node embeddings for clustering  

% Apply k-means clustering to the embeddings
[cluster_indexs, centroids] = kmeans(embeddings, k);

end