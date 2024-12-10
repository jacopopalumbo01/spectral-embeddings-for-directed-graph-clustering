function [eigval_red,eigvec_red] = BCS(W, P, k, plot, verbose, graph_name)
% BCS - Perform cyclic spectral clustering
%
%% Syntax:
%        BCS(W, P, k, plot, verbose, graph_name)
%
%% Input Arguments:
%       *Required Input Arguments*
%       - W:                        Adjacency matrix
%       - P:                        Transition matrix
%       - k:                        Number of clusters
%       - plot:                     Plot the eigenvalues (true/false)
%       - verbose:                  Verbose output 
%
%% Output:
%       - eigval_red:               Cycle eigenvalues
%       - eigvec_red:               Cycle eigenvectors
%


if nargin < 4
    plot = true;
    verbose = true;
end

if nargin < 5
    verbose = true;
end


%% Find the k cycle eigenvalues
% Perform Arnoldi iteration
%[Q, H] = Arnoldi(P, 10);
  
% Remove last row from H
%H(size(H,1),:) = [];

% Compute Eigenvectors and Eigenvalues
%[V, D] = eig(H);
[V, D] = eig(P);

% Get first k Eigenvalues and Eigenvectors
D = diag(D);
modulus = abs(D);

% Select the k largest cycle eigenvalues based on modulus
eigval_red = zeros(k,1);
eigvec_red = zeros(size(W,1),k);
for i = 1:k
    [M, j] = max(modulus);
    eigval_red(i) = D(j);
    eigvec_red(:,i)  = V(:,j);

    V(:,j) = [];
    D(j) = [];
    modulus(j) = [];
end

if verbose
    fprintf('Eigvec matrix with rows: %d and cols: %d\n',...
        size(eigvec_red,1),size(eigvec_red,2));
end

if plot
    % Plot eigenvalues and eigenvectors 
    if nargin < 6
        PlotCyclicEig(D,eigval_red,eigvec_red);
    else
        PlotCyclicEig(D,eigval_red,eigvec_red, graph_name);
    end
end


end