clear;close all;

%% Parameters
seeds = [0; 1; 2; 5; 42; 123; 1234; 1991; 2001; 2025];

ks= [2; 3; 4; 5; 6; 7; 8; 9; 10];   % Number of clusters
n = 5000;                           % Number of nodes
p = 0.008;                          % Connection prob. inside cluster
q = p;                              % Connection prob. inter-clusters
mu = 0.00;                          % Noise parameter

for nc = 1:size(ks,1)
    k = ks(nc);
    fprintf("----------\n")
    fprintf("| k = %d |\n", k);
    fprintf("----------\n")
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
    
    for i = 1:size(seeds,1)
        fprintf("Seed: %d\n", seeds(i));
        rng(seeds(i));
        fprintf("   Generating DSBM...\n");
        [W, labels] = generate_dag_dsbm(k, p, q, c, F);
        fprintf("   Saving graph...\n");
        save(sprintf("synthetic/generated/DSBM_%dblocks_%dnodes_%fnoise_%dseed.mat", k, sum(c,"all"), mu, seeds(i)), "W", "labels");
    end
end
