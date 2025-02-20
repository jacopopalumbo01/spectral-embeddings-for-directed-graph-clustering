clear all; clc;

addpaths_SC;
warning off;
rng(42);
plot_graph_flag = false;
plot_eig_flag   = false;
verbose         = false; 

% Synthetically generate an unperturbed uniform block-cyclic
graph_name = "Unperturbed Uniform Block-Acyclic graph";
fprintf("%s\n", graph_name);
n   = 100;          % Number of nodes
k   = 8;            % Number of blocks  
P   = zeros(k,k);   % Block connection probability

% Block membership distribution
rho = ones(k,1);
rho = rho / sum(rho);

% Set block connection probability
conn_prob = 0.5;
for i = 1:k
    for j = 1:k
        if i < j
            P(i,j) = conn_prob;
        end
    end
end

%% Generation of unperturbed graph
fprintf("   Generating unperturbed graph\n")
[W, nodes] = StochasticBlockmodel(n,k,rho,P);

%% Estimate number of blocks
fprintf("   Estimating number of blocks\n");
k_estimate = EstimateNumBlocksAcyclic(W,n/10);
fprintf("   Estimated: %d blocks. Ground truth: %d blocks\n", k_estimate, k);

% Plot the graph
if plot_graph_flag == 1
    PlotCyclic(W, k, transpose(nodes), graph_name)
end

% Get clusters
[cluster_index, ~] = BAS(W, k, plot_eig_flag, verbose, graph_name);


%% Evaluation
fprintf("----------\nEvaluation\n----------\n")
[RCut, NCut, NMI, FScore] = ComputeMetrics(nodes,cluster_index,W);

fprintf("   RCut:       %f\n", RCut);
fprintf("   NCut:       %f\n", NCut);
fprintf("   NMI:        %f\n", NMI);
fprintf("   F-score:    %f\n", FScore);

fprintf("\n\n");


% Synthetically generate an unperturbed not uniform block-cycle
graph_name = "Unperturbed not Uniform Block-Acyclic graph";
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
        if i < j
            P(i,j) = conn_prob;
        end
    end
end

%% Generation of unperturbed graph
fprintf("   Generating unperturbed graph\n")
[W, nodes] = StochasticBlockmodel(n,k,rho,P);

%% Estimate number of blocks
fprintf("   Estimating number of blocks\n");
k_estimate = EstimateNumBlocksAcyclic(W,n/10);
fprintf("   Estimated: %d blocks. Ground truth: %d blocks\n", k_estimate, k);

% Plot the input graph
if plot_graph_flag == 1    
    PlotCyclic(W, k, transpose(nodes), graph_name)
end
% Get clusters
[cluster_index, ~] = BAS(W, k, plot_eig_flag, verbose, graph_name);

%% Evaluation
fprintf("----------\nEvaluation\n----------\n")
[RCut, NCut, NMI, FScore] = ComputeMetrics(nodes,cluster_index,W);

fprintf("   RCut:       %f\n", RCut);
fprintf("   NCut:       %f\n", NCut);
fprintf("   NMI:        %f\n", NMI);
fprintf("   F-score:    %f\n", FScore);

fprintf("\n\n");


% Synthetically generate a perturbed block-acyclic graph
graph_name = "Perturbed Uniform Block-Acyclic graph";
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
        if i < j
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

%% Estimate number of blocks
fprintf("   Estimating number of blocks\n");
k_estimate = EstimateNumBlocksAcyclic(W,n/10);
fprintf("   Estimated: %d blocks. Ground truth: %d blocks\n", k_estimate, k);

% Plot the graph
if plot_graph_flag == 1
    PlotCyclic(W, k, transpose(nodes), graph_name);
end

% Get clusters
[cluster_index, ~] = BAS(W, k, plot_eig_flag, verbose, graph_name);

%% Evaluation
fprintf("----------\nEvaluation\n----------\n")
[RCut, NCut, NMI, FScore] = ComputeMetrics(nodes,cluster_index,W);

fprintf("   RCut:       %f\n", RCut);
fprintf("   NCut:       %f\n", NCut);
fprintf("   NMI:        %f\n", NMI);
fprintf("   F-score:    %f\n", FScore);

fprintf("\n\n");


% Synthetically generate a perturbed block-acyclic
graph_name = "Perturbed not Uniform Block-Acyclic graph";
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
        if i < j
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

%% Estimate number of blocks
fprintf("   Estimating number of blocks\n");
k_estimate = EstimateNumBlocksAcyclic(W,n/10);
fprintf("   Estimated: %d blocks. Ground truth: %d blocks\n", k_estimate, k);

% Plot the graph
if plot_graph_flag == 1
    PlotCyclic(W, k, transpose(nodes), graph_name);
end

% Get clusters
[cluster_index, ~] = BAS(W, k, plot_eig_flag, verbose, graph_name);

%% Evaluation
fprintf("----------\nEvaluation\n----------\n")
[RCut, NCut, NMI, FScore] = ComputeMetrics(nodes,cluster_index,W);

fprintf("   RCut:       %f\n", RCut);
fprintf("   NCut:       %f\n", NCut);
fprintf("   NMI:        %f\n", NMI);
fprintf("   F-score:    %f\n", FScore);

