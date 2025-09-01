clear; close all;

k= 25;                               % Number of clusters
n = 2500;                           % Number of nodes
p = 0.008;                          % Connection prob. inside cluster
q = p;                              % Connection prob. inter-clusters
mu = 0.00;                          % Noise parameter

% Base value for each cluster
base_val = floor(n / k);
c = base_val * ones(k, 1);

% Distribute the remainder
remainder = n - sum(c);
c(1:remainder) = c(1:remainder) + 1;

% Build Cluster level orientation probabilities for Directed Acyclic DSBM
F = zeros(k,k);
for i = 1:k
    for j = 1:k
        if j == i + 1 || j == i + 2
            F(i,j) = mu;
        elseif j == i - 1 || j == i - 2
            F(i,j) = 1 - mu;
        else
            F(i,j) = 0.5;
        end
    end
end

% Generate graph
[W, ~] = generate_dag_dsbm(k, p, q, c, F);

% Plot sparsity
spy(W, 'k.', 10);
axis off;