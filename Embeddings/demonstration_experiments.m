function ratio = compute_distance_ratio(X)
    % X: Input data matrix (n x d), n points in d dimensions
    n = size(X, 1);
    
    % Compute pairwise Euclidean distances
    D = pdist2(X, X);  % n x n distance matrix
    
    % Set diagonal to Inf for nearest neighbor calculation
    D_temp = D;
    D_temp(1:n+1:end) = Inf;
    
    % Compute average distance to nearest neighbor (across all points)
    nn_dists = min(D_temp, [], 2);  % Nearest neighbor distances
    avg_nn = mean(nn_dists);
    
    % Compute average distance to random point (global mean)
    D_temp2 = D;
    D_temp2(1:n+1:end) = 0;        % Set diagonal to 0
    total_dist = sum(D_temp2(:));   % Sum of all pairwise distances
    num_pairs = n*(n-1)/2;          % Number of unique pairs
    avg_random = total_dist / (2*num_pairs);  % Adjust for double-counting
    
    ratio = avg_random / avg_nn;
end

%{
% Demo: Test with low vs. high dimensions
n = 500;  % Number of points

% Low-dimensional case (d=1)
d_low = 1;
X_low = randn(n, d_low);
R_low = compute_distance_ratio(X_low);

% High-dimensional case (d=1000)
d_high = 1000;
X_high = randn(n, d_high);
R_high = compute_distance_ratio(X_high);

fprintf('d = %d: R = %.4f (nearest neighbors are much closer)\n', d_low, R_low);
fprintf('d = %d: R = %.4f (distances concentrated → curse of dimensionality)\n', d_high, R_high);

% Visualization: How R changes with dimension
dims = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000];
ratios = zeros(size(dims));

for i = 1:length(dims)
    X = randn(n, dims(i));
    ratios(i) = compute_distance_ratio(X);
end

figure;
semilogx(dims, ratios, 'bo-', 'LineWidth', 1.5);
hold on;
yline(1, 'r--', 'LineWidth', 1.5);
xlabel('Dimension (d)');
ylabel('Ratio R');
title('Curse of Dimensionality: R vs. Dimension');
legend('Actual R', 'R=1 (Theoretical Limit)', 'Location', 'southeast');
grid on;
%}

clear;

clusters = [2;3;4;5;6;7;8;9;10];
ratio = zeros(size(clusters,1),1);
ratio_tSNE = zeros(size(clusters,1),1);

for i = 1:size(clusters,1)
    % Load adjacency
    loaded = load(sprintf("DSBM_%dblocks_5000nodes_0.000000noise.mat", clusters(i)));
    W      = loaded.W;
    labels = loaded.labels;
            
    % Get number of labels
    k = size(unique(labels), 1);
    
    [embeddings, ~] = SVD_unscaled_embeddings(W,k);
    ratio(i) = compute_distance_ratio(embeddings);

    embed_tsne = tsne(embeddings);
    ratio_tSNE(i) = compute_distance_ratio(embed_tsne);
end

figure;
loglog(clusters, ratio, 'bo-', 'LineWidth', 1.5);
hold on;
loglog(clusters, ratio_tSNE, 'go-', 'LineWidth', 1.5);
hold on;
yline(1, 'r--', 'LineWidth', 1.5);
xlabel('Clusters (k)');
ylabel('Ratio R');
title('Curse of Dimensionality: R vs. No. Clusters k');
legend('Actual R', 'R after tSNE', 'R=1 (Theoretical Limit)', 'Location', 'southeast');