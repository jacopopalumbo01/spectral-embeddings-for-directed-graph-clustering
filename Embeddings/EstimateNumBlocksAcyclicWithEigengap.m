function [k, P, V, Eigvals] = EstimateNumBlocksAcyclicWithEigengap(W, max_nclust)
% EstimateNumBlocksAcyclicWithEigengap - Compute number of clusters given 
% the adjacency matrix using eigengaps criterion.
%
% Input:
%   - W:            Adjacency matrix
%   - max_nclust:   Maximum number of clusters to consider
%
% Output:
%   - k:            Estimated number of clusters
%   - P:            Transition probability matrix
%   - V:            Eigenvectors of P
%   - D:            Eigenvalues of P

% Compute the transition matrix
P       = AcyclicTransitionMatrix(W);
% Small epsilon to avoid division by zero
eps_val = 1e-10; 
% Compute eigenvalues and eigenvectors
% ARPACK, Implicitly Restarted Arnoldi Method (IRAM), O(k*n) per iteration.
% [V, D]   = eigs(P,max_nclust,'lm'); 
% LAPACK DGEEV/ZGEEV, QR algorithm with Hessenberg reduction, total O(n^3).
if ~issparse(W)
    [V, D]   = eig(full(P)); 
else
    [V, D]   = eigs(P, 2*max_nclust);
end

Eigvals = diag(D);

% Apply filtering criteria to eigvals
% valid_indices = (real(eigvals) < 1) & (imag(eigvals) >= 0);
valid_indices = (imag(Eigvals) >= -eps_val);
Eigvals = Eigvals(valid_indices);

% Sort eigenvalues by descending absolute value, and
% accordingly the eignvectors
[~, sort_idx] = sort(abs(Eigvals), 'descend');
Eigvals = Eigvals(sort_idx);
V       = V(:,sort_idx); 

% Keep only the first k_to_find eigenvalues
k_to_find = min(max_nclust, length(Eigvals));
k_to_find = max(2, k_to_find); % Ensure at least 2 eigenvalues

% Truncate eigenvalues
Eigvals = Eigvals(1:k_to_find);

% Compute all relative eigengaps

eiggaps = zeros(length(Eigvals)-1, 1);

for i = 2:length(Eigvals)-1
    eiggaps(i-1) = (abs(Eigvals(i)) - abs(Eigvals(i+1))) / (abs(Eigvals(i+1)) + eps_val);
end

% Select candidate k with maximum eiggap score
[~, optimal_k_idx] = max(abs(eiggaps));
% +1 because eiggaps start from the 2nd eigenvalue
k = optimal_k_idx + 1;  

% Plot eigenvalues and eigengaps for inspection
figure;

subplot(1,2,1);
plot(1:length(Eigvals), abs(Eigvals), 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Index');
ylabel('|\lambda|');
title('Sorted Eigenvalues');
grid on;

subplot(1,2,2);
plot(2:length(Eigvals), abs(eiggaps), 'ro-', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
plot(optimal_k_idx+1, abs(eiggaps(optimal_k_idx)), 'ks', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Index');
ylabel('Relative Eigengap');
title('Relative Eigengap Spectrum');
legend('Eigengap', 'Selected k');
grid on;

end
