clear;close all;
rng(1991);

%% Parameters
k = 5;                  % Number of clusters
p = 0.008;              % Connection prob. inside cluster
q = p;                  % Connection prob. inter-clusters
c = ones(k,1) * 1000;   % Nodes per cluster
mu = 0.0;                 % Noise parameter

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

[W, labels] = generate_dag_dsbm(k, p, q, c, F);

save(sprintf("DSBM_%dblocks_%dnodes", k, sum(c,"all")), "W", "labels");

%{
% Visualize
G = digraph(W);
figure;
plot(G);
title('Directed Acyclic DSBM Graph');


%}