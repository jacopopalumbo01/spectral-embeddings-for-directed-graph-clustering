function [modularity] = Modularity(W, membership)
% Modularity - Returns the modularity associated with the provided directed
%              graph and clusters.
%
% Input:
%   - W:            The adjacency matrix.
%   - membership:   The cluster membership function (vector)
%
% Output:
%   - modularity:   The modularity score

% Number of clusters
k = max(membership);

% Total number of edges
num_edges = sum(W, 'all');

% Compute the in-degree and out-degree per cluster
k_out = zeros(k, 1);
k_in = zeros(k, 1);

for i = 1:k
    k_out(i) = sum(W(membership == i, :), 'all'); % Sum of outgoing edges from cluster i
    k_in(i)  = sum(W(:, membership == i), 'all'); % Sum of incoming edges to cluster i
end

% Compute the actual fraction of edges within clusters
num_edges_i2j = zeros(k, k);
for i = 1:k
    for j = 1:k
        num_edges_i2j(i,j) = sum(W(membership == i, membership == j), 'all') / num_edges;
    end
end

% Compute expected edge distribution
expected_edges = (k_out * k_in') / (num_edges^2); % Expected edge distribution

% Compute modularity using only diagonal elements
modularity = sum(diag(num_edges_i2j) - diag(expected_edges)); 

end
