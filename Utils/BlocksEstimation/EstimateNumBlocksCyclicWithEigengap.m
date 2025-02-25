function [k] = EstimateNumBlocksCyclicEigengap(W, max_nclust)
% EstimateNumBlocksCyclicEigengap - Compute number of clusters given 
% the adjacency matrix using eigengaps criterion.
%
% Input:
%   - A:            Adjacency matrix
%   - max_k         Maximum number of clusters to consider
%
% Output:
%   - k:            Estimated number of clusters

%% Step 1: Compute the transition matrix
P = TransitionMatrix(W);

%% Step 2: Compute eigenvalues and eigenvectors
[~, D] = eig(P); 
eigvals = diag(D);

% Apply the filtering criteria
valid_indices = (real(eigvals) < 1) & (imag(eigvals) >= 0);
eigvals = eigvals(valid_indices);

% Sort eigenvalues by descending absolute value so that the largest (in modulus) come first
[~, sort_idx] = sort(abs(eigvals), 'descend');
eigvals = eigvals(sort_idx);

% Keep only the first k_to_find eigenvalues/eigenvectors.
% Here we define k_to_find as the smaller of floor(max_k/2) and the available number.
k_to_find = min(floor(max_nclust/2), length(eigvals));
% (Ensure at least 2 for a valid embedding.)
k_to_find = max(2, k_to_find);

% Truncate the eigenpairs
eigvals = eigvals(1:k_to_find);


%% Step 3: Compute all the eigengaps
eiggaps = zeros(size(eigvals,1)-1, 1);

for i = 2:size(eigvals,1)
    eiggaps(i-1) = (eigvals(i - 1) - eigvals(i)) / eigvals(i);
end

%% Step 4: Select candidate k with maximum eiggap score
[~, optimal_k_idx] = max(abs(eiggaps));
k = optimal_k_idx + 1;  % because candidates start at 2

end
