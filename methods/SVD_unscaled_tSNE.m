function [unscaled_adj, clusters] = SVD_unscaled_tSNE(A, k)
% Code for Clustering Using SVD as per: https://arxiv.org/pdf/1108.2228
% Input: A adjacency matrix
%        k num_blocks

if ~issparse(A)
    % Step 1: Compute the Singular Value Decomposition
    [U, S, V] = svd(A);
    
    % Step 2: Reduce Dimensionality
    % Select the top-k singular vectors
    U_k = U(:, 1:k);   % Left singular vectors
    V_k = V(:, 1:k);   % Right singular vectors
    S_k = S(1:k, 1:k); % Top-k singular values
else
    [U_k, S_k, V_k] = svds(A);
end

% Create unscaled adjacency spectral embedding
unscaled_adj = [U_k, V_k];

% Reduce dimensionality
embeddings = tsne(unscaled_adj);

% Apply k-means clustering to the embeddings
[clusters, ~] = kmeans(embeddings, k,'Replicates', 20);

end