function A = incidence_to_adjacency_unsymm(Inc)
% Convert an incidence matrix to a sparse adjacency matrix
% Inc - incidence matrix [n_edges x n_nodes]

% Get the number of nodes and edges
[n_edges, n_nodes] = size(Inc);

% Initialize row, column, and value arrays to store non-zero entries
rows   = zeros(n_edges,1);
cols   = zeros(n_edges,1);
values = zeros(n_edges,1);

% Loop through each edge (row in the incidence matrix)
for i = 1:n_edges
    % Find the non-zero entries in the incidence matrix row
    nodes_from = find(Inc(i, :) == -1);
    nodes_to   = find(Inc(i, :) == 1);

    % If there are exactly two nodes involved in the edge
    % For directed graphs, set both A(from,to) to 1
    rows(i)   = nodes_from;
    cols(i)   = nodes_to;
    values(i) = 1;
end

% Create the sparse adjacency matrix using row, column, and value arrays
A = sparse(rows, cols, values, n_nodes, n_nodes);
end
