clear all; clc;

addpaths_SC;
rng(42);
% Synthetically generate an unperturbed uniform block-cyclic
graph_name = "Unperturbed Uniform Block-Cyclic graph";
fprintf("%s\n", graph_name);
n   = 100;          % Number of nodes
k   = 8;            % Number of blocks  
P   = zeros(k,k);   % Block connection probability

% Block membership distribution
rho = ones(k,1);
rho = rho / sum(rho);

% Set block connection probability
conn_prob = 0.7;
for i = 1:k
    for j = 1:k
        if i + 1 == j || (i == k && j == 1)
            P(i,j) = conn_prob;
        end
    end
end

%% Generation of unperturbed graph
fprintf("   Generating unperturbed graph\n")
[W, nodes] = StochasticBlockmodel(n,k,rho,P);


% Plot the graph
PlotCyclic(W, k, nodes, graph_name)

% Get clusters
[cluster_index, ~] = BCS(W, k, true, true, graph_name);


%% Evaluation
fprintf("----------\nEvaluation\n----------\n")
[RCut, NCut, NMI, FScore] = ComputeMetrics(nodes,cluster_index,W);

fprintf("   RCut:       %f\n", RCut);
fprintf("   NCut:       %f\n", NCut);
fprintf("   NMI:        %f\n", NMI);
fprintf("   F-score:    %f\n", FScore);

fprintf("\n\n");


% Synthetically generate an unperturbed uniform block-cycle
graph_name = "Unperturbed not Uniform Block-Cycle graph";
fprintf("%s\n", graph_name);
n   = 100;          % Number of nodes
k   = 8;            % Number of blocks  
P   = zeros(k,k);   % Block connection probability

% Block membership distribution
rho = [0.18; 0.2; 0.05; 0.2; 0.14; 0.04; 0.07; 0.13];

% Set block connection probability
conn_prob = 0.7;
for i = 1:k
    for j = 1:k
        if i + 1 == j || (i == k && j == 1)
            P(i,j) = conn_prob;
        end
    end
end

%% Generation of unperturbed graph
fprintf("   Generating unperturbed graph\n")
[W, nodes] = StochasticBlockmodel(n,k,rho,P);


% Plot the graph
PlotCyclic(W, k, nodes, graph_name)

% Get clusters
[cluster_index, ~] = BCS(W, k, true, true, graph_name);

%% Evaluation
fprintf("----------\nEvaluation\n----------\n")
[RCut, NCut, NMI, FScore] = ComputeMetrics(nodes,cluster_index,W);

fprintf("   RCut:       %f\n", RCut);
fprintf("   NCut:       %f\n", NCut);
fprintf("   NMI:        %f\n", NMI);
fprintf("   F-score:    %f\n", FScore);

fprintf("\n\n");


% Synthetically generate a perturbed block-cyclic graph
graph_name = "Perturbed Uniform Block-Cyclic graph";
fprintf("%s\n", graph_name);
n   = 100;          % Number of nodes
k   = 8;            % Number of blocks  
P   = zeros(k,k);   % Block connection probability
epsilon = 0.4;      % Perturbation magnitude

% Block membership distribution
rho = ones(k,1);
rho = rho / sum(rho);

% Set block connection probability
conn_prob = 0.7;
for i = 1:k
    for j = 1:k
        if i + 1 == j || (i == k && j == 1)
            P(i,j) = conn_prob;
        end
    end
end

%% Generation of unperturbed graph
fprintf("   Generating unperturbed graph\n")
[W, nodes] = StochasticBlockmodel(n,k,rho,P);

%% Generation of perturbing graph
% Create perturbing block connection probability matrix
Q  = rand(k,k);

% Computing new perturbing block connection probability matrix given
% the magnitude of perturbation epsilon
P2 = epsilon * Q;

% Generate perturbing graph
fprintf("   Generating perturbing graph\n");
[A, ~] = StochasticBlockmodel(n,k,rho,P2);

%% Combine unperturbed and perturbing graph
fprintf("   Combining unperturbed and perturbing graph\n");
W = CombineBlockmodels(W,A);

% Plot the graph
PlotCyclic(W, k, nodes, graph_name);

% Get clusters
[cluster_index, ~] = BCS(W, k, true, true, graph_name);

%% Evaluation
fprintf("----------\nEvaluation\n----------\n")
[RCut, NCut, NMI, FScore] = ComputeMetrics(nodes,cluster_index,W);

fprintf("   RCut:       %f\n", RCut);
fprintf("   NCut:       %f\n", NCut);
fprintf("   NMI:        %f\n", NMI);
fprintf("   F-score:    %f\n", FScore);

fprintf("\n\n");


% Synthetically generate a perturbed block-cycle
graph_name = "Perturbed not Uniform Block-Cycle graph";
fprintf("%s\n", graph_name);
n   = 100;          % Number of nodes
k   = 8;            % Number of blocks  
P   = zeros(k,k);   % Block connection probability
epsilon = 0.4;      % Perturbation magnitude

% Block membership distribution
rho = [0.18; 0.2; 0.05; 0.2; 0.14; 0.04; 0.07; 0.13];

% Set block connection probability
conn_prob = 0.7;
for i = 1:k
    for j = 1:k
        if i + 1 == j || (i == k && j == 1)
            P(i,j) = conn_prob;
        end
    end
end

%% Generation of unperturbed graph
fprintf("   Generating unperturbed graph\n")
[W, nodes] = StochasticBlockmodel(n,k,rho,P);

%% Generation of perturbing graph
% Create perturbing block connection probability matrix
Q  = rand(k,k);

% Computing new perturbing block connection probability matrix given
% the magnitude of perturbation epsilon
P2 = epsilon * Q;

% Generate perturbing graph
fprintf("   Generating perturbing graph\n");
[A, ~] = StochasticBlockmodel(n,k,rho,P2);

%% Combine unperturbed and perturbing graph
fprintf("   Combining unperturbed and perturbing graph\n");
W = CombineBlockmodels(W,A);

% Plot the graph
PlotCyclic(W, k, nodes, graph_name);

% Get clusters
[cluster_index, ~] = BCS(W, k, true, true, graph_name);

%% Evaluation
fprintf("----------\nEvaluation\n----------\n")
[RCut, NCut, NMI, FScore] = ComputeMetrics(nodes,cluster_index,W);

fprintf("   RCut:       %f\n", RCut);
fprintf("   NCut:       %f\n", NCut);
fprintf("   NMI:        %f\n", NMI);
fprintf("   F-score:    %f\n", FScore);

