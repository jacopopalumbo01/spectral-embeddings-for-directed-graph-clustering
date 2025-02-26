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
    P = TransitionMatrix(W);
    
    % Compute eigenvalues and eigenvectors
    [V, EigVals] = eigs(P, k, 'lm');
    EigVals      = diag(EigVals); % Extract eigenvalues
    
    if verbose
        fprintf('Selected %d largest eigenvalues and their eigenvectors.\n', k);
    end
    
    % Plot eigenvalues and eigenvectors
    if plotFlag
        D = eig(P);
        modulus = abs(D);
            
        % Select k largest eigenvalues and corresponding eigenvectors
        % Note: We are interested in the eigenvalues with largest modulus
        [~, indices] = maxk(modulus, k);
        cycle_eigvals = D(indices);
        D(indices) = [];
        
        PlotCyclicEig(D, cycle_eigvals, V, "");
    end
    
    % Combine real and imaginary parts of eigenvectors
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
        [cluster_indexs, centroids] = kmeans(data_normalized, k, 'Distance', 'sqeuclidean','Replicates', 20);
    else
        % Perform k-means clustering on real and imaginary parts
        [cluster_indexs, centroids] = kmeans(data_real_imag, k, 'Distance', 'sqeuclidean','Replicates', 20);
    end
    
    if verbose
        fprintf('Clustering completed. %d clusters identified.\n', k);
    end
end
    