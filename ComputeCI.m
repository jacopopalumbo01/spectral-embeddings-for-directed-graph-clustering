function [CI_vals, CIsz_vals, CIvol_vals, TopCIvol, TopCIsz, TopTF] = ComputeCI(W, clusters, c)
% ComputeCI
% ----------
% Computes the pairwise Cut Imbalance (CI), size-weighted (CI^sz),
% and volume-weighted (CI^vol) metrics between all cluster pairs.
%
% Also computes:
%   - TopCIvol: Sum of top c CI^vol scores.
%   - TopCIsz : Sum of top c CI^sz scores.
%   - TopTF   : Sum of top c Trade Flow scores.
%
% INPUTS:
%   W         : Adjacency matrix (n x n), directed.
%   clusters  : Cluster assignment vector (n x 1), values in {1,...,k}.
%   c         : Number of top pairs to select (e.g. c = 2*(k-1) for acyclic, c = k for cylcic).
%
% OUTPUTS:
%   CI_vals   : Array of unscaled CI values for all cluster pairs.
%   CIsz_vals : CI scaled by min cardinality of the cluster pair.
%   CIvol_vals: CI scaled by min volume of the cluster pair.
%   TopCIvol  : Sum of top-c largest CI^vol values.
%   TopCIsz   : Sum of top-c largest CI^sz values.
%   TopTF     : Sum of top-c largest Trade Flow values.

k = max(clusters);

CI_vals    = [];
CIsz_vals  = [];
CIvol_vals = [];
TF_vals    = [];

% Loop through all unique cluster pairs (i < j)
for i = 1:k
    Ai = find(clusters == i);  % Nodes in cluster i
    for j = i+1:k
        Aj = find(clusters == j);  % Nodes in cluster j

        % Directed edge weights between clusters
        w_ij = sum(W(Ai, Aj), 'all');  % Edges from i -> j
        w_ji = sum(W(Aj, Ai), 'all');  % Edges from j -> i
        w_total = w_ij + w_ji;         % Total weight across the cut

        % CI value (cut imbalance) \in [0, 0.5]
        CI = 0.5 * abs(w_ij - w_ji) / max(w_total, eps);

        % Cluster sizes (cardinality)
        size_i = length(Ai);
        size_j = length(Aj);

        % Cluster volumes (in + out degrees)
        vol_i = sum(W(Ai,:), 'all') + sum(W(:,Ai), 'all');
        vol_j = sum(W(Aj,:), 'all') + sum(W(:,Aj), 'all');

        % Scaled CI metrics
        CIsz  = CI * min(size_i, size_j);
        CIvol = CI * min(vol_i, vol_j);

        % Trade Flow metric (directional imbalance)
        TF = abs(w_ij - w_ji);

        % Append to lists
        CI_vals(end+1)    = CI;
        CIsz_vals(end+1)  = CIsz;
        CIvol_vals(end+1) = CIvol;
        TF_vals(end+1)    = TF;
    end
end

% Sort descending to get top-c contributions
CIvol_sorted = sort(CIvol_vals, 'descend');
CIsz_sorted  = sort(CIsz_vals,  'descend');
TF_sorted    = sort(TF_vals,    'descend');

% Sum the top-c values (clip if fewer pairs exist)
c = min(c, length(CIvol_sorted));
TopCIvol = sum(CIvol_sorted(1:c));
TopCIsz  = sum(CIsz_sorted(1:c));
TopTF    = sum(TF_sorted(1:c));

end
