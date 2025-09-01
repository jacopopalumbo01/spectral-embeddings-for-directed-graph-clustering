%{
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


clusters = [2;3;4;5;6;7;8;9;10];
ratio = zeros(size(clusters,1),1);
ratio_tSNE = zeros(size(clusters,1),1);

for i = 1:size(clusters,1)
    % Load adjacency
    loaded = load(sprintf("DSBM_%dblocks_5000nodes_0.000000noise_0seed.mat", clusters(i)));
    W      = loaded.W;
    labels = loaded.labels;
            
    % Get number of labels
    k = size(unique(labels), 1);
    
    [embeddings, ~] = SVD_unscaled(W,k);
    ratio(i) = compute_distance_ratio(embeddings);

    embed_tsne = tsne(embeddings);
    ratio_tSNE(i) = compute_distance_ratio(embed_tsne);
end
%}
figure;
semilogx(clusters, ratio, 'bx-', 'LineWidth', 4, 'MarkerSize', 20);
hold on;
yline(1, 'r--', 'LineWidth', 4);
xlabel('Clusters (k)', "Interpreter","latex");
ylabel('ADR', "Interpreter","latex");
legend('Actual ADR', 'ADR=1 (Theoretical Limit)', "Interpreter", "Latex");
legend("boxoff");
set(gca, "fontsize", 50);

figure;
loglog(clusters, ratio, 'bx-', 'LineWidth', 4, 'MarkerSize', 20);
hold on;
loglog(clusters, ratio_tSNE, 'gx-', 'LineWidth', 4, 'MarkerSize', 20);
hold on;
yline(1, 'r--', 'LineWidth', 4);
xlabel('Clusters (k)', "Interpreter", "Latex");
ylabel('ADR', "Interpreter", "Latex");
legend('Actual ADR', 'ADR after tSNE', 'ADR=1 (Theoretical Limit)', "Interpreter", "Latex");
legend("boxoff");
set(gca, "fontsize", 50);