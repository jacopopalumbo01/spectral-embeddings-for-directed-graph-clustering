function [cluster_indeces, centroids] = BASnew(W, k, transition, beta, plotFlag, verbose, graph_name)
% BASnew - Perform acyclic spectral clustering
%
% Syntax:
%        [cluster_indexs, centroids] = BAS(W, k, transition, beta, plotFlag, verbose, graph_name)
%
% Input Arguments:
%       - W (required):            Adjacency matrix (NxN)
%       - k (required):            Number of clusters
%       - transition (required):   Type of transition matrix
%                                  ("transition"/"power")
%       - beta (required):         Parameter used when transition=power
%       - plotFlag (optional):     Plot eigenvalues (true/false, default: false)
%       - verbose (optional):      Verbose output (true/false, default: false)
%       - graph_name (optional):   Name of the graph (default: '')
%
% Output:
%       - cluster_indeces:          Clustering index labels
%       - centroids:               Centroids obtained during clustering

% Validate inputs
if nargin < 3, transition = "transition"; end
if nargin < 5, plotFlag = false; end
if nargin < 6, verbose = false; end
if nargin < 7, graph_name = ''; end


if ~ismatrix(W) || size(W, 1) ~= size(W, 2)
    error('Input W must be a square adjacency matrix.');
end
if k <= 0 || k > size(W, 1)
    error('Input k must be a positive integer and <= number of nodes.');
end
 
% Compute transition probability matrix
if transition == "transition"
    P = AcyclicTransitionMatrix(W);
elseif transition == "power"
    P = AdjustableAcyclicTransitionMatrix(W,beta);
end


% Compute eigenvalues and eigenvectors
[V, EigVals] = eigs(P,k,'lm', "maxit", 10000);
EigVals = diag(EigVals); % Extract eigenvalues

% Filter out eigenvalues
valid_indices = (real(EigVals) < 1) & (imag(EigVals) >= 0);
EigVals = EigVals(valid_indices);
V = V(:, valid_indices);

k_to_find = min(floor(k/2));
EigVals = EigVals(1:k_to_find);
V = V(:, 1:k_to_find);

if verbose
    fprintf('Selected %d largest eigenvalues and their eigenvectors.\n', k);
end

% Plot eigenvalues and eigenvectors.
% Note that if we are computing only the k largest ones,
% the first argument should be empty []
if plotFlag
    D = eig(P);
    modulus = abs(D);
        
    % Select k largest eigenvalues and corresponding eigenvectors
    % Note: We are interested in the eigenvalues with largest modulus
    [~, indices] = maxk(modulus, k);
    cycle_eigvals = D(indices);
    D(indices) = [];
    
    PlotCyclicEig(D, cycle_eigvals, V, "");
    Eigs_All = eig(P);
end

% Combine real and imaginary parts of eigenvectors
%data_real_imag = [real(cycle_eigvecs), imag(cycle_eigvecs)];
data_real_imag = [real(V), imag(V)];

norm_eigvecs = 0;
if norm_eigvecs == 1
    % Step 1: Compute norms and detect zero rows
    norms = vecnorm(data_real_imag, 2, 2);  % Row-wise norms
    zero_rows = norms == 0;                % Identify zero rows
    norms(zero_rows) = 1;                  % Avoid division by zero

    % Step 2: Normalize non-zero rows
    data_normalized = data_real_imag ./ norms;

    % Step 3: Handle zero rows 
    data_normalized(zero_rows, :) = 0;
    % Perform k-means clustering on real and imaginary parts
    [cluster_indeces, centroids] = kmeans(data_normalized, k, 'Distance', 'sqeuclidean','Replicates', 20);
else
    % Perform k-means clustering on real and imaginary parts
    [cluster_indeces, centroids] = kmeans(data_real_imag, k, 'Distance', 'sqeuclidean','Replicates', 20);
end

if verbose
    fprintf('Clustering completed. %d clusters identified.\n', k);
end
end
    