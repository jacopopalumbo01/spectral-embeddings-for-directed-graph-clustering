function score = PrecedenceAlignmentScore(W, clusters)
% Edge Direction Consistency (Precedence Alignment Score)
% The fraction of edges that respect the cluster order

    [u, v] = find(W);
    total_edges = length(u);
    consistent = sum(clusters(u) <= clusters(v)); % Respect order
    score = consistent / total_edges;
end
