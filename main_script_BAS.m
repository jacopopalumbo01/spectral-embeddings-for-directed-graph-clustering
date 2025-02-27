clear all; close all;

addpaths_SC;
warning off;
rng(42);

%% Flags
plot_graph_flag = false;
plot_eig_flag   = false;
verbose         = false; 

%% Graph properties
n         = 100;          % Number of nodes
k         = 8;            % Number of blocks  
conn_prob = 0.7;          % Connection probability between blocks

% Uniform block membership distribution
rho_uniform = ones(k,1);
rho_uniform = rho_uniform / sum(rho_uniform);

% Not uniform block membership distribution
rho_not_uniform = [0.18; 0.2; 0.05; 0.2; 0.14; 0.04; 0.07; 0.13];

epsilon = 0.4;            % Perturbation magnitude

% Generate block connection probability
P = ConnectionProbabilityMatrix("acyclic", k, conn_prob);

results = zeros(4,5); % Where we store all the results

%% Synthetically generate an unperturbed uniform block-cyclic
graph_name = "Unperturbed Uniform Block-Acyclic graph";
fprintf("%s\n", graph_name);

% Generate block connection probability
P = ConnectionProbabilityMatrix("acyclic", k, conn_prob);

% Generation of unperturbed graph
fprintf("   Generating unperturbed graph\n")
[W, nodes] = StochasticBlockmodel(n,k,rho_uniform,P);

% Estimate number of blocks
fprintf("   Estimating number of blocks\n");
k_estimate = EstimateNumBlocksAcyclic(W,n/10);
fprintf("   Estimated: %d blocks. Ground truth: %d blocks\n", k_estimate, k);

% Plot the graph
if plot_graph_flag == 1
    PlotCyclic(W, k, transpose(nodes), graph_name)
end

% Get clusters
[cluster_index, ~] = BAS(W, k, "transition", 1, plot_eig_flag, verbose, graph_name);


% Evaluation
[RCut, NCut, NMI, FScore, modularity] = ComputeMetrics(nodes,cluster_index,W);
results(1,:) = [RCut, NCut, NMI, FScore, modularity];


%% Synthetically generate an unperturbed not uniform block-cycle
graph_name = "Unperturbed not Uniform Block-Acyclic graph";
fprintf("%s\n", graph_name);

% Generation of unperturbed graph
fprintf("   Generating unperturbed graph\n")
[W, nodes] = StochasticBlockmodel(n,k,rho_not_uniform,P);

% Estimate number of blocks
fprintf("   Estimating number of blocks\n");
k_estimate = EstimateNumBlocksAcyclic(W,n/10);
fprintf("   Estimated: %d blocks. Ground truth: %d blocks\n", k_estimate, k);

% Plot the input graph
if plot_graph_flag == 1    
    PlotCyclic(W, k, transpose(nodes), graph_name)
end
% Get clusters
[cluster_index, ~] = BAS(W, k, "transition", 1, plot_eig_flag, verbose, graph_name);

% Evaluation
[RCut, NCut, NMI, FScore, modularity] = ComputeMetrics(nodes,cluster_index,W);
results(2,:) = [RCut, NCut, NMI, FScore, modularity];


%% Synthetically generate a perturbed block-acyclic graph
graph_name = "Perturbed Uniform Block-Acyclic graph";
fprintf("%s\n", graph_name);

% Generate block connection probability
P = ConnectionProbabilityMatrix("acyclic", k, conn_prob);

% Generation of unperturbed graph
fprintf("   Generating unperturbed graph\n")
[W, nodes] = StochasticBlockmodel(n,k,rho_uniform,P);

% Generation of perturbing graph
% Create perturbing block connection probability matrix
Q  = rand(k,k);

% Computing new perturbing block connection probability matrix given
% the magnitude of perturbation epsilon
P2 = epsilon * Q;

% Generate perturbing graph
fprintf("   Generating perturbing graph\n");
[A, ~] = StochasticBlockmodel(n,k,rho_uniform,P2);

% Combine unperturbed and perturbing graph
fprintf("   Combining unperturbed and perturbing graph\n");
W = CombineBlockmodels(W,A);

% Estimate number of blocks
fprintf("   Estimating number of blocks\n");
k_estimate = EstimateNumBlocksAcyclic(W,n/10);
fprintf("   Estimated: %d blocks. Ground truth: %d blocks\n", k_estimate, k);

% Plot the graph
if plot_graph_flag == 1
    PlotCyclic(W, k, transpose(nodes), graph_name);
end

% Get clusters
[cluster_index, ~] = BAS(W, k, "transition", 1, plot_eig_flag, verbose, graph_name);

% Evaluation
[RCut, NCut, NMI, FScore, modularity] = ComputeMetrics(nodes,cluster_index,W);
results(3,:) = [RCut, NCut, NMI, FScore, modularity];


%% Synthetically generate a perturbed block-acyclic
graph_name = "Perturbed not Uniform Block-Acyclic graph";
fprintf("%s\n", graph_name);

% Generate block connection probability
P = ConnectionProbabilityMatrix("acyclic", k, conn_prob);

% Generation of unperturbed graph
fprintf("   Generating unperturbed graph\n")
[W, nodes] = StochasticBlockmodel(n,k,rho_not_uniform,P);

% Generation of perturbing graph
% Create perturbing block connection probability matrix
Q  = rand(k,k);

% Computing new perturbing block connection probability matrix given
% the magnitude of perturbation epsilon
P2 = epsilon * Q;

% Generate perturbing graph
fprintf("   Generating perturbing graph\n");
[A, ~] = StochasticBlockmodel(n,k,rho_not_uniform,P2);

% Combine unperturbed and perturbing graph
fprintf("   Combining unperturbed and perturbing graph\n");
W = CombineBlockmodels(W,A);

% Estimate number of blocks
fprintf("   Estimating number of blocks\n");
k_estimate = EstimateNumBlocksAcyclic(W,n/10);
fprintf("   Estimated: %d blocks. Ground truth: %d blocks\n", k_estimate, k);

% Plot the graph
if plot_graph_flag == 1
    PlotCyclic(W, k, transpose(nodes), graph_name);
end

% Get clusters
[cluster_index, ~] = BAS(W, k, "transition", 1, plot_eig_flag, verbose, graph_name);

% Evaluation
[RCut, NCut, NMI, FScore, modularity] = ComputeMetrics(nodes,cluster_index,W);
results(4,:) = [RCut, NCut, NMI, FScore, modularity];

%% Report the results
fprintf("----------------------\n");
fprintf("@  Graph Properties  @\n")
fprintf("----------------------\n");

fprintf("Nodes:                     %d\n", n);
fprintf("Blocks:                    %d\n", k);
fprintf("Connection probability:    %f\n", conn_prob);
fprintf("Perturbation Magnitude:    %f\n\n", epsilon);

fprintf("----------------------------\n");
fprintf("@          Results         @\n")
fprintf("----------------------------\n\n");

T = array2table(results, "VariableNames", ["RCut", "NCut", "NMI", "F-Score", "Modularity"], ...
    "RowNames",["Unperturbed Uniform Block-Cyclic", "Unperturbed Not Uniform Block-Cycle", ...
    "Perturbed Uniform Block-Cyclic", "Perturbed Not Uniform Block-Cycle"]);
disp(T);