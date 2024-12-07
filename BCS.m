function [eigval_red,eigvec_red] = BCS(W, P, k, plot, pedantic)
% Perform cyclic spectral clustering
% Input
% W: Adjacency matrix
% P: Transition matrix
% k: Number of clusters
% plot: Plot the eigenvalues (true/false)
% pedantic: Pedantic output (true/false)
% Output
% [eigenvalues, eigenvectors]


if nargin < 4
    plot = true;
    pedantic = true;
end

if nargin < 5
    pedantic = true;
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

if pedantic
    fprintf('Eigvec matrix with rows: %d and cols: %d\n',...
        size(eigvec_red,1),size(eigvec_red,2));
end

if plot
    % Plot eigenvalues and eigenvectors 
    PlotCyclicEig(D,eigval_red,eigvec_red);
end


end