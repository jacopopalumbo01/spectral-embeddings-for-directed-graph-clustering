function [embeddings, clusters] = BCS_embeddings(W, k)
% BCS_embeddings - Perform acyclic spectral clustering
%
% Syntax:
%        [embeddings, centroids] = BCS_embeddings(W, k)
%
% Input Arguments:
%       - W (required):            Adjacency matrix (NxN)
%       - k (required):            Number of clusters
%
% Output:
%       - embeddings:              Obtained embeddings
%       - clusters:                Inferred labels
 

% Compute transition probability matrix
P = TransitionMatrix(W);


% Compute eigenvalues and eigenvectors
[V, EigVals] = eig(P);
EigVals = diag(EigVals); % Extract eigenvalues

% Sort by largest magnitude                           
[~, idx] = sort(abs(EigVals), 'descend');     
EigVals = EigVals(idx);               
V = V(:, idx);  

% Filter out eigenvalues
valid_indices = (real(EigVals) <= 1) & (imag(EigVals) >= 0);

EigVals = EigVals(valid_indices);
V = V(:, valid_indices);

k_to_find = min(floor(k/2));

EigVals = EigVals(1:k_to_find);
V = V(:, 1:k_to_find);

% Combine real and imaginary parts of eigenvectors
embeddings = [real(V), imag(V)];

% Apply k-means clustering to the embeddings
[clusters, ~] = kmeans(embeddings, k,'Replicates', 20);

end
    