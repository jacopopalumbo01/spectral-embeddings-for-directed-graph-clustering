function A = incidence_to_adjacency_unsymm_large(Inc)
    % Get the number of edges (rows) and nodes (columns)
    [n_edges, n_nodes] = size(Inc);

    % Find all source (-1) and destination (+1) indices
    [row_idx_src, src_nodes] = find(Inc == -1); % Edge rows and their sources
    [row_idx_dst, dst_nodes] = find(Inc == 1);  % Edge rows and their destinations

    % Ensure correct edge matching by row index
    [sorted_rows_src, sort_idx_src] = sort(row_idx_src);
    [sorted_rows_dst, sort_idx_dst] = sort(row_idx_dst);

    % Apply the same sorted order to src_nodes and dst_nodes
    src_nodes = src_nodes(sort_idx_src);
    dst_nodes = dst_nodes(sort_idx_dst);

    % Debugging: Print a sample of correctly matched edges
    disp('Sample edges from Inc (first 10):');
    disp([src_nodes(1:10), dst_nodes(1:10)]);

    % Construct adjacency matrix preserving multiple edges
    A = sparse(src_nodes, dst_nodes, ones(size(src_nodes)), n_nodes, n_nodes);

    % Debugging Output
    disp(['nnz(A) (after fix): ', num2str(nnz(A))]);
end
