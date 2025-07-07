function avg_conductance = Compute_Conductance(W, clusters)
% Compute average conductance over all clusters in a partition
% For directed graphs: computes outgoing volume and cut edges
%
% Inputs:
%   - W        : (n x n) adjacency matrix
%   - clusters : (n x 1) cluster assignments
%
% Output:
%   - avg_conductance: mean conductance across clusters

    k = max(clusters);
    n = size(W,1);
    conductances = zeros(1, k);

    for c = 1:k
        S = find(clusters == c);
        Sc = find(clusters ~= c);

        % Outgoing cut edges from S to not-S
        cut = sum(W(S, Sc), 'all');

        % Volume = total outgoing edges from S
        vol = sum(W(S, :), 'all');

        % Conductance: cut(S, ~S) / vol(S)
        conductances(c) = cut / max(vol, eps);
    end

    avg_conductance = mean(conductances);
end
