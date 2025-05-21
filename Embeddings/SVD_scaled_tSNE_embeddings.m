function [scaled_adj, clusters] = SVD_scaled_tSNE_embeddings(A, k)
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

% Create scaled adjacency spectral embedding
scaled_adj = [U_k * sqrt(S_k), V_k * sqrt(S_k)];

% Reduce dimensionality
embeddings = tsne(scaled_adj);

% Apply k-means clustering to the embeddings
[clusters, ~] = kmeans(embeddings, k,'Replicates', 20);

end