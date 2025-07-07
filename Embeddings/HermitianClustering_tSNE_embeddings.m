function [G, clusters] = HermitianClustering_tSNE_embeddings(W, k)
% HermitianClustering - Spectral clustering for directed graphs
% using Hermitian adjacency matrix as described in the paper.
%
% Input:
%   W : (n x n) adjacency matrix (binary or weighted)
%   k : desired number of clusters
%
% Output:
%   clusters : (n x 1) cluster assignments

    n = size(W, 1);
    
    % Step 1: Construct Hermitian matrix H
    H = ConstructHermitianMatrix(W);

    % Step 2: Determine l (number of eigenvectors)
    if mod(k, 2) == 0
        l = k;
    else
        l = k - 1;
    end

    % Step 3: Compute the l eigenvectors with largest magnitude eigenvalues
    [V, D] = eig(H, 'vector');
    [~, idx] = sort(abs(D), 'descend');
    V_selected = V(:, idx(1:l));

    % Step 4: Compute projection matrix P = sum_j v_j * v_j^*
    P = zeros(n, n);
    for j = 1:l
        P = P + V_selected(:, j) * V_selected(:, j)';
    end

    % Step 5: Second eigen decomposition of P to get embedding G
    [G, ~] = eigs(P, l, 'largestreal', 'Tolerance', 1e-6);

    G = real(G);

    % Reduce dimensionality
    embeddings = tsne(G);

    % Step 6: Run k-means on rows of G
    clusters = kmeans(embeddings, k, 'Replicates', 10, 'MaxIter', 300);
end
