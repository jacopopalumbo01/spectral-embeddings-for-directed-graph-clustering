function [cluster_indexs,centroids] = BCS(W, k, plot, verbose, graph_name)
% BCS - Perform cyclic spectral clustering
%
%% Syntax:
%        BCS(W, P, k, plot, verbose, graph_name)
%
%% Input Arguments:
%       *Required Input Arguments*
%       - W:                        Adjacency matrix
%       - k:                        Number of clusters
%       - plot:                     Plot the eigenvalues (true/false)
%       - verbose:                  Verbose output 
%
%% Output:
%       - cluster_index:           Clustering index labels
%       - centroids:               Centroids obtained during clustering
%


if nargin < 3
    plot = true;
    verbose = true;
end

if nargin < 4
    verbose = true;
end

%% Compute the transition probability matrix
P = TransitionMatrix(W);

if plot
    % Plot sparsity pattern of W
    if nargin < 5
        figure;
        spy(W,'k.',15)
        axis off
    else
        figure;
        spy(W,'k.',15)
        axis off
        title(sprintf("%s Adjacency Matrix", graph_name));
    end
    
end

%% Find the k cycle eigenvalues

% Compute Eigenvectors and Eigenvalues
[V, D] = eig(P);

% Get first k Eigenvalues and Eigenvectors
D = diag(D);
modulus = abs(D);

% Select the k largest cycle eigenvalues based on modulus
cycle_eigvals = zeros(k,1);
cycle_eigvecs = zeros(size(W,1),k);
for i = 1:k
    [M, j] = max(modulus);
    cycle_eigvals(i) = D(j);
    cycle_eigvecs(:,i)  = V(:,j);

    V(:,j) = [];
    D(j) = [];
    modulus(j) = [];
end

if verbose
    fprintf('Eigvec matrix with rows: %d and cols: %d\n',...
        size(cycle_eigvecs,1),size(cycle_eigvecs,2));
end

if plot
    % Plot eigenvalues and eigenvectors 
    if nargin < 5
        PlotCyclicEig(D,cycle_eigvals,cycle_eigvecs);
    else
        PlotCyclicEig(D,cycle_eigvals,cycle_eigvecs, graph_name);
    end
end


% Extract the real and imaginary part 
% from the cycle eigenvectors
cycle_real = real(cycle_eigvecs);
cycle_imag = imag(cycle_eigvecs);
% The new data matrix is [num_nodes, 2xecycle_eigenvecs]
data_real_imag = [cycle_real, cycle_imag];

%% Step 3: K-means on the rows
[cluster_indexs, centroids] = kmeans(data_real_imag, k, 'Distance', 'sqeuclidean');

end