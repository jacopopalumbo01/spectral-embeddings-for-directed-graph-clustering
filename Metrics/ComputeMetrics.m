function [RCut, NCut, NMI, FScore, modularity] = ComputeMetrics(true_membership, inferred_membership, W)
% ComputeMetrics - Compute all the metrics used during evaluation
%
%% Syntax:
%        ComputeMetrics(clusters, nodes)
%
%% Input Arguments:
%       *Required Input Arguments*
%       - true_membership:          True membership of nodes
%       - inferred_membership:      Inferred membership of nodes
%       - W:                        Adjacency matrix
%
%% Output:
%       - RCut:                     Ratio Cut
%       - NCut:                     Normalized Cut
%       - NMI:                      Normalized Mutual Information
%       - FScore:                   F-Score
%


%% Internal metrics (based only on graph structure and weights)
normalized = 1;
NCut = computeRCutValue(inferred_membership,W,normalized);
RCut = computeRCutValue(inferred_membership,W,~normalized);

%% External metrics (based on labels)
% Compute nmi
NMI = nmi(true_membership, inferred_membership);

[inferred_labels,~] = label_data(inferred_membership,true_membership,1);

% Compute f-score
[Scores] = evaluate_scores(true_membership,inferred_labels);
FScore = Scores(3);

modularity = Compute_modularity(W, inferred_membership);
end