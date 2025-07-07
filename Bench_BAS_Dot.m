%% Main Benchmark Script for BAS on Directed Graphs in DOT Format
clear; clc; close all;
rng(1991);

addpaths_SC(); % Make sure your BAS functions are in path

%% Parameters
graph_cases = {
% 'instance_CG_N5_K4_nzP0d4_random_comm_weights.dot';
% 'instance_CG_N6_K5_nzP0d5_random_comm_weights.dot';
% 'instance_CG_N7_K2_nzP0d6_random_comm_weights.dot';
% 'instance_CG_N7_K7_nzP0d2_random_comm_weights.dot';
% 'instance_CG_N8_K3_nzP0d5_random_comm_weights.dot';
% 'instance_CG_N9_K5_nzP0d2_random_comm_weights.dot';
% 'instance_exp_N10_K8_nzP0d2_random_comm_weights.dot';
% 'instance_exp_N15_K4_nzP0d2_random_comm_weights.dot';
% 'instance_kNN_N10_K8_nzP0d2_random_comm_weights.dot';
% 'instance_kNN_N15_K4_nzP0d25_random_comm_weights.dot';
% 'instance_pregel_cc_gyro_m_random_comm_weights.dot';
% 'instance_simple_pagerank_gyro_m_random_comm_weights.dot';
% 'instance_snni_graphchallenge_1024neurons_120layers_random_comm_weights.dot';
% 'instance_spmv_N25_nzP0d2_random_comm_weights.dot';
% 'instance_spmv_N35_nzP0d18_random_comm_weights.dot';
% 'instance_spmv_N40_nzP0d15_random_comm_weights.dot';
'RandomBand_p80_b5_100_419_A_random_comm_weights.dot';
};

weight_type = 'none'; % Options: 'comm_weight', 'work_weight', 'mem_weight', 'none'

num_trials       = 1;             % How many times to run BAS on each graph
transition_type  = 'transition';  % 'transition' or 'power'
beta             = 0.1;           % Only if transition == 'power'
plotFlag         = false;         % Plot eigenvalues inside BAS
verbose          = false;         % BAS verbosity
estimate_blocks  = true;          % Whether to estimate k or fix it manually

%% Initialize Results
nc = length(graph_cases);
results = struct();

fprintf('========================================================\n');
fprintf('      Block-Acyclic Spectral Clustering (BAS) Benchmark \n');
fprintf('========================================================\n\n');

fprintf('Graph files input    : DOT format\n');
fprintf('Selected weight type : %s\n', weight_type);
fprintf('Transition type      : %s\n', transition_type);
if strcmp(transition_type, 'power')
    fprintf('Beta parameter       : %.4f\n', beta);
end
fprintf('Trials per graph     : %d\n\n', num_trials);

%% Loop over cases
for c = 1:nc
    fprintf('--------------------------------------------------------\n');
    fprintf('Processing graph file: %s\n', graph_cases{c});
    
    %% Step 1: Load graph from .dot file
    [W, node_weights] = DOTtoAdjacencyMatrix(graph_cases{c}, weight_type);
    
    n_nodes = size(W, 1);
    n_edges = nnz(W);
    
    fprintf('Nodes   : %d\n', n_nodes);
    fprintf('Edges   : %d\n', n_edges);
    
    %% Initialize metrics storage for this graph
    RCut_all       = zeros(num_trials, 1);
    NCut_all       = zeros(num_trials, 1);
    Modularity_all = zeros(num_trials, 1);
    k_est_all      = zeros(num_trials, 1);
    time_all       = zeros(num_trials, 1);
    
    %% Step 2: Run BAS clustering multiple times
    for t = 1:num_trials
        fprintf('Trial %d/%d...\n', t, num_trials);
        tStart = tic;
        
        % Estimate k or set manually
        if estimate_blocks
            k_est = EstimateNumBlocksAcyclicWithEigengap(W, floor(n_nodes / 10));
        else
            error('You must specify k if estimate_blocks is false.');
        end
        
        % Run BAS clustering
        [clusters, centroids] = BAS(W, k_est, transition_type, beta, plotFlag, verbose, graph_cases{c});
        
        %% Compute partitioning metrics
        normalized = 1;
        NCut       = computeRCutValue(clusters, W, normalized);
        RCut       = computeRCutValue(clusters, W, ~normalized);
        Modularity = Compute_modularity(W, clusters);
        
        %% Store metrics
        RCut_all(t)       = RCut;
        NCut_all(t)       = NCut;
        Modularity_all(t) = Modularity;
        k_est_all(t)      = k_est;
        time_all(t)       = toc(tStart);
    end
    
    %% Step 3: Save results for this graph
    results(c).graph_name  = graph_cases{c};
    results(c).n_nodes     = n_nodes;
    results(c).n_edges     = n_edges;
    results(c).RCut        = RCut_all;
    results(c).NCut        = NCut_all;
    results(c).Modularity  = Modularity_all;
    results(c).k_est       = k_est_all;
    results(c).time        = time_all;
    
    %% Step 4: Print per-graph summary
    fprintf('Results for graph: %s\n', graph_cases{c});
    fprintf('--------------------------------------------------------\n');
    fprintf('Nodes        : %d\n', n_nodes);
    fprintf('Edges        : %d\n', n_edges);
    fprintf('Estimated k  : %.2f (mean across trials)\n', mean(k_est_all));
    fprintf('RCut (avg)   : %.4f | Std: %.4f\n', mean(RCut_all), std(RCut_all));
    fprintf('NCut (avg)   : %.4f | Std: %.4f\n', mean(NCut_all), std(NCut_all));
    fprintf('Modularity   : %.4f | Std: %.4f\n', mean(Modularity_all), std(Modularity_all));
    fprintf('Time (avg)   : %.4fs | Std: %.4fs\n', mean(time_all), std(time_all));
    fprintf('--------------------------------------------------------\n\n');
    
    %% Step 5: Visualization for each graph
    % Plot clusters (pass clusters as the node assignment, consistent with PlotCyclic)
    PlotCyclic(W, k_est, clusters);
end

%% Step 6: Print all results summary
fprintf('======================= Summary ========================\n');
fprintf('%-35s %-6s %-6s %-8s %-10s %-10s %-10s\n', ...
    'Graph', 'Nodes', 'Edges', 'Blocks', 'NCut', 'Modularity', 'Time(s)');

for c = 1:nc
    fprintf('%-35s %-6d %-6d %-8.2f %-10.4f %-10.4f %-10.4f\n', ...
        results(c).graph_name, ...
        results(c).n_nodes, ...
        results(c).n_edges, ...
        mean(results(c).k_est), ...
        mean(results(c).NCut), ...
        mean(results(c).Modularity), ...
        mean(results(c).time));
end
fprintf('========================================================\n');
