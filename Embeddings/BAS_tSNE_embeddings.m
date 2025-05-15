function [embeddings, clusters] = BAS_embeddings(W, k)
% BAS_tSNE_embeddings - Perform acyclic spectral clustering with
% dimensionality reduction
%
% Syntax:
%        [embeddings, clusters] = BAS_tSNE_embeddings(W, k)
%
% Input Arguments:
%       - W (required):            Adjacency matrix (NxN)
%       - k (required):            Number of clusters
%
% Output:
%       - embeddings:              Obtained embeddings
%       - clusters:                Inferred clusters
 
% Compute transition probability matrix
P = AcyclicTransitionMatrix(W);


% Compute eigenvalues and eigenvectors
[V, EigVals] = eigs(P,k,'lm', "maxit", 10000);
EigVals = diag(EigVals); % Extract eigenvalues

% Filter out eigenvalues
valid_indices = (real(EigVals) <= 1) & (imag(EigVals) >= 0);

EigVals = EigVals(valid_indices);
V = V(:, valid_indices);

k_to_find = min(floor(k/2));
EigVals = EigVals(1:k_to_find);
V = V(:, 1:k_to_find);

% Combine real and imaginary parts of eigenvectors
embeddings = [real(V), imag(V)];

% Reduce dimensionality
red_embed = tsne(embeddings);

% Apply k-means clustering to the embeddings
[clusters, ~] = kmeans(red_embed, k,'Replicates', 20);

end
    