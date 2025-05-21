function [W, labels] = generate_dag_dsbm(k, p, q, c, F)
% generate_dag_dsbm - Generates DSBM as in 10.1109/BigData55660.2022.10020413
%
% Syntax:
%        [W, labels] = generate_dag_dsbm(k, p, q, c, F)
%
% Input Arguments:
%       - k (required):             Number of clusters.
%       - p (required):             Probability that two vertices in the 
%                                   same cluster have an edge between them.
%       - q (required):             Probability that two vertices in different
%                                   clusters have an edge between them.
%       - c (required):             A vector of length k whose entries are 
%                                   the number of vertices in each cluster.
%       - F (required):             Cluster level orientation probabilities,
%                                   size(F) must be (k,k).
%
% Output:
%       - W:                        Adjacency matrix
%       - labels:                   Nodes labels
 
% Validate input
if size(c,1) ~= k
    error('Length of c must be equal to k.');
end
if ~isequal(size(F), [k, k])
    error('F must be a k x k matrix.');
end

n = sum(c);              % total number of nodes
labels = zeros(n,1);     % initialize labels
W = zeros(n);            % initialize adjacency matrix

% Compute indices for each cluster
cluster_idx = cell(k,1);

start_idx = 1;
for i = 1:k
    cluster_idx{i} = start_idx:(start_idx + c(i) - 1);
    labels(cluster_idx{i}) = i;
    start_idx = start_idx + c(i);
end

% Fill adjacency matrix
for i = 1:k
    for j = 1:k
        idx_i = cluster_idx{i};
        idx_j = cluster_idx{j};

        n_i = length(idx_i);
        n_j = length(idx_j);

        if i == j
            % Intra-cluster edges
            A = rand(n_i);  % generate a dense matrix

            % Determine where to put an edge
            mask = A < p;

            % Remove all self-loops (set diagonal = 0)
            diagonalMask = logical(eye(size(mask)));
            mask(diagonalMask) = 0;
            
            % Save the edges
            W(idx_i, idx_i) = double(mask);
        else
            % Inter-cluster edges
            prob = q * F(i,j);
            mask = rand(n_i, n_j) < prob;
            W(idx_i, idx_j) = double(mask);
        end
    end
end
end
