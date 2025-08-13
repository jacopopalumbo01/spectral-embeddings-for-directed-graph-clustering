function [NMI, FScore] = Compute_ext_metrics(true_membership, inferred_membership)
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
%       - NMI:                      Normalized Mutual Information
%       - FScore:                   F-Score
%

%% External metrics (based on labels)
[inferred_labels,~] = label_data(inferred_membership,true_membership,1);

% Compute nmi
NMI = nmi(true_membership, inferred_labels);

% Compute f-score
[Scores] = evaluate_scores(true_membership,inferred_labels);
FScore = Scores(3);
end