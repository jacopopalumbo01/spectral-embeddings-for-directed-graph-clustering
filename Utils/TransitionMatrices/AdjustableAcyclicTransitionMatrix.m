function P = AdjustableAcyclicTransitionMatrix(W, beta)
% AdjustableAcyclicTransitionMatrix - Compute the transition matrix from an 
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

% Prepare rows with out-degree == 0. We want 1 in each column of these rows
% to perform 1/size(W,1) when normalizing.
W(row_sums == 0, :) = 1;

% Recompute row sums after modification
row_sums = sum(W, 2);

% Power factor
row_sums = row_sums .^ beta;

% Normalize each row to create the transition matrix
P = W ./ row_sums;

end
