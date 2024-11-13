function [tau] = BCS(W, k)

    %% Step 1: Compute the transition matrix P
    P = TransitionMatrix(W);

    %% Step 2: Find the k cycle eigenvalues
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

    for i = 1:k
        [M, j] = max(modulus);
        eigval(i, :) = D(j,:);
        eigvec(:,i) = V(:,j);
        V(:,j) = [];
        modulus(j) = [];
    end
    size(eigvec)
    %eigval
    %eigvec

    %% Step 3: K-means on the rows
    [tau, C] = LloydCluster(eigvec, k, 100)
    %[tau, C] = kmeans(eigvec, k);
end