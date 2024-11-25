function [tau] = BCS(W, P, k)

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
    fprintf('Eigvec matrix with rows: %d and cols: %d\n',...
        size(eigvec_red,1),size(eigvec_red,2));

    % Plot eigenvalues
    PlotCyclicEig(D,eigval_red,eigvec_red);

    %% Step 3: K-means on the rows
    %[tau, C] = LloydCluster(eigvec_red, k, 100);
    [tau, C] = LloydClusterSensitive(eigvec_red, k, 500);
end