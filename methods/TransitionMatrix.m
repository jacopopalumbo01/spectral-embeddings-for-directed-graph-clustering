function P = TransitionMatrix(W)
% TransitionMatrix - Compute the transition matrix from an adjacency matrix
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

% Prevent division by zero by setting zero rows to uniform distribution
row_sums(row_sums == 0) = 1;

% Normalize each row to create the transition matrix
P = W ./ row_sums;

% Set rows with zero out-degree explicitly to zeros (optional)
zero_rows = (sum(W, 2) == 0);
P(zero_rows, :) = 0;

end
