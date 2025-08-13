function P = AcyclicTransitionMatrix(W)
% TransitionMatrix - Compute the transition matrix from an adjacency matrix
% According to formula
% P_{ij} =
% \begin{cases}
% \frac{1}{d_i^{\text{out}}} W_{ij} & \text{if } d_i^{\text{out}} > 0, \\
% \frac{1}{n} & \text{otherwise}.
% \end{cases}
%
% Syntax:
%        P = TransitionMatrix(W)
%
% Input:
%        W - Adjacency matrix (NxN)
%
% Output:
%        P - Transition matrix (NxN)

% Compute row sums (out-degrees)
row_sums = sum(W, 2);

% Prepare rows with out-degree == 0. We want 1 in each column of these rows
% to perform 1/size(W,1) when normalizing.
W(row_sums == 0, :) = 1;

% Recompute row sums after modification
row_sums = sum(W, 2);

% Normalize each row to create the transition matrix
P = W ./ row_sums;

end
