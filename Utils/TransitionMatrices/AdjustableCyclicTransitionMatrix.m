function P = AdjustableCyclicTransitionMatrix(W, beta)
% AdjustableCyclicTransitionMatrix - Compute the transition matrix from an 
% adjacency matrix given a power factor
%
% Input:
%   - W:        Adjacency matrix
%   - beta:     Power factor
%
% Output:
%   - P:        Transition matrix

% Compute row sums (out-degrees)
row_sums = sum(W, 2);

% Prevent division by zero by setting zero rows to uniform distribution
row_sums(row_sums == 0) = 1;

% Power factor
row_sums = row_sums .^ beta;

% Normalize each row to create the transition matrix
P = W ./ row_sums;

% Set rows with zero out-degree explicitly to zeros (optional)
zero_rows = (sum(W, 2) == 0);
P(zero_rows, :) = 0;

end
