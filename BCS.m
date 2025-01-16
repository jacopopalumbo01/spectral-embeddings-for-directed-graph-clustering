function [cluster_indexs, centroids] = BCS(W, k, plotFlag, verbose, graph_name)
    % BCS - Perform cyclic spectral clustering
    %
    % Syntax:
    %        [cluster_indexs, centroids] = BCS(W, k, plotFlag, verbose, graph_name)
    %
    % Input Arguments:
    %       - W (required):            Adjacency matrix (NxN)
    %       - k (required):            Number of clusters
    %       - plotFlag (optional):     Plot eigenvalues (true/false, default: false)
    %       - verbose (optional):      Verbose output (true/false, default: false)
    %       - graph_name (optional):   Name of the graph (default: '')
    %
    % Output:
    %       - cluster_indexs:          Clustering index labels
    %       - centroids:               Centroids obtained during clustering
    
    % Validate inputs
    if nargin < 3, plotFlag = false; end
    if nargin < 4, verbose = false; end
    if nargin < 5, graph_name = ''; end
    
    if ~ismatrix(W) || size(W, 1) ~= size(W, 2)
        error('Input W must be a square adjacency matrix.');
    end
    if k <= 0 || k > size(W, 1)
        error('Input k must be a positive integer and <= number of nodes.');
    end
    
    % Compute transition probability matrix
    P = TransitionMatrix_dms(W);
    
    % Plot adjacency matrix sparsity pattern
    if plotFlag
        figure;
        spy(W, 'k.', 15);
        axis off;
        if ~isempty(graph_name)
            title(sprintf("%s Adjacency Matrix", graph_name));
        end
    end
    
    % Compute eigenvalues and eigenvectors
    [V, D] = eig(P);
    D = diag(D); % Extract eigenvalues
    modulus = abs(D);
    
    % Select k largest eigenvalues and corresponding eigenvectors
    % Note: We are interested in the eigenvalues with largest modulus
    [~, indices] = maxk(modulus, k);
    cycle_eigvals = D(indices);
    cycle_eigvecs = V(:, indices);
    
    if verbose
        fprintf('Selected %d largest eigenvalues and their eigenvectors.\n', k);
    end
    
    % Plot eigenvalues and eigenvectors
    if plotFlag
        PlotCyclicEig(D, cycle_eigvals, cycle_eigvecs, graph_name);
    end
    
    % Combine real and imaginary parts of eigenvectors
    data_real_imag = [real(cycle_eigvecs), imag(cycle_eigvecs)];
    
    % Perform k-means clustering
    disp('Here in kmeans');
    [cluster_indexs, centroids] = kmeans(data_real_imag, k, 'Distance', 'sqeuclidean','Replicates', 20);
    
    if verbose
        fprintf('Clustering completed. %d clusters identified.\n', k);
    end
end
    