%% Main benchmark script for synthetic graphs
clear; close all;
rng(1991);

%addpaths_DirAC();

%% Parameters
graph_cases = {
    % n = 100, decay = 0.7
    % '2blocks_100nodes.mat';
    % '3blocks_100nodes.mat';
    % '4blocks_100nodes.mat';
    % '5blocks_100nodes.mat';
    % '6blocks_100nodes.mat';
    % '7blocks_100nodes.mat';
    % '8blocks_100nodes.mat';
    % '9blocks_100nodes.mat';
    % '10blocks_100nodes.mat';
    %
    %  % n = 200, decay = 0.7
    % '2blocks_200nodes.mat';
    % '3blocks_200nodes.mat';
    % '4blocks_200nodes.mat';
    % '5blocks_200nodes.mat';
    % '6blocks_200nodes.mat';
    % '7blocks_200nodes.mat';
     '8blocks_200nodes.mat';
     '8blocks_200nodes_unsorted.mat';
    % '9blocks_200nodes.mat';
    % '10blocks_200nodes.mat';
    %
    % % n = 500, decay = 0.6
    % '2blocks_500nodes.mat';
    % '3blocks_500nodes.mat';
    % '4blocks_500nodes.mat';
    % '5blocks_500nodes.mat';
    % '6blocks_500nodes.mat';
    % '7blocks_500nodes.mat';
    % '8blocks_500nodes.mat';
    % '9blocks_500nodes.mat';
    % '10blocks_500nodes.mat';
    %
    % % n = 1000, decay = 0.5
    % '2blocks_1000nodes.mat';
    % '3blocks_1000nodes.mat';
    % '4blocks_1000nodes.mat';
    % '5blocks_1000nodes.mat';
    % '6blocks_1000nodes.mat';
    % '7blocks_1000nodes.mat';
    % '8blocks_1000nodes.mat';
    % '9blocks_1000nodes.mat';
    % '10blocks_1000nodes.mat';

    % 'test_mat_perm.mat';

    % '8-blocks_2500-nodes_0.000000-pert.mat';
    % '8-blocks_2500-nodes_0.500000-pert.mat';
    % '8-blocks_2500-nodes_0.900000-pert.mat';
    % 'eco.mat'
    };

weight_type      = 'none';
trans_reduce     = 0;
power_fac        = 2;
num_trials       = 1;
transition_type  = 'transition';
beta             = 0.1;
plotFlag         = 0;
verbose          = false;
estimate_blocks  = 0;

%% Initialize Results
nc = length(graph_cases);
results = struct();

fprintf('========================================================\n');
fprintf('      Block-Acyclic Spectral Clustering (BAS) Benchmark \n');
fprintf('========================================================\n\n');

if strcmp(transition_type, 'power')
    fprintf('Beta parameter       : %.4f\n', beta);
end
fprintf('Trials per graph     : %d\n\n', num_trials);

%% Loop over cases
for c = 1:nc
    fprintf('--------------------------------------------------------\n');
    fprintf('Processing graph file: %s\n', graph_cases{c});

    [~, ~, ext] = fileparts(graph_cases{c});

    switch lower(ext)
        case '.dot'
            [W_orig, node_weights] = DOTtoAdjacencyMatrix(graph_cases{c}, weight_type);
            W = W_orig;
        case '.mat'
            mat_data = load(graph_cases{c});
            W = mat_data.W;
            labels = mat_data.labels;
            W_orig = W;
        otherwise
            error('Unsupported file extension: %s', ext);
    end

    if trans_reduce == 1
        G = digraph(W);
        H = transreduction(G);
        W = adjacency(H);
    end

    n_nodes = size(W_orig, 1);
    n_edges = nnz(W_orig);
    fprintf('Nodes   : %d\n', n_nodes);
    fprintf('Edges   : %d\n', n_edges);

    % Initialize result storage for each method
    methods = {'bas', 'svd', 'skew', 'herm'};
    % methods = {'svd'};
    for m = 1:length(methods)
        method = methods{m};
        TopCIvol_all.(method) = zeros(num_trials, 1);
        TopCIsz_all.(method)  = zeros(num_trials, 1);
        TopTF_all.(method)    = zeros(num_trials, 1);
        Phi_all.(method)      = zeros(num_trials, 1);
        NMI_all.(method)      = zeros(num_trials, 1);
        PaS_all.(method)      = zeros(num_trials, 1);
        Fscore_all.(method)   = zeros(num_trials, 1);
        time_all.(method)     = zeros(num_trials, 1);
    end
    k_est_all = zeros(num_trials, 1);

    k = max(unique(labels));
    for t = 1:num_trials

        k_est_all(t) = k;
        c_thres = 2*(k-1);

        % % Clustering methods
        t_start = tic;
        eig_method    = 2;
        project_space = 0;

        %clusterings.dirac = DiRAC(W, k, eig_method, power_fac, project_space);
        %time_all.dirac(t) = toc(t_start);


        t_start = tic;
        clusterings.bas = BCSnew(W, k);
        time_all.bas(t) = toc(t_start);

        t_start = tic;
        clusterings.svd   = SVD_clustering_suss(full(W), k);
        time_all.svd(t) = toc(t_start);

        t_start = tic;
        clusterings.skew  = SkewSymmetricClustering(W, k);
        time_all.skew(t) = toc(t_start);

        t_start = tic;
        clusterings.herm  = HermitianClustering(W, k);
        time_all.herm(t) = toc(t_start);

        for m = 1:length(methods)
            method = methods{m};
            tStart = tic;
            clusters = clusterings.(method);

            Phi = Compute_Conductance(W, clusters);
            PaS = PrecedenceAlignmentScore(W, clusters);
            [CI_vals, CIsz_vals, CIvol_vals, TopCIvol, TopCIsz, TopTF] = ComputeCI(W, clusters, c_thres);
            [NMI, Fscore] = Compute_ext_metrics(labels, clusters);            

            TopCIvol_all.(method)(t) = TopCIvol;
            TopCIsz_all.(method)(t)  = TopCIsz;
            TopTF_all.(method)(t)    = TopTF;
            PaS_all.(method)(t)      = PaS;
            Phi_all.(method)(t)      = Phi;
            NMI_all.(method)(t)      = NMI;
            Fscore_all.(method)(t)   = Fscore;
        end
    end

    %% Save per-graph result
    for m = 1:length(methods)
        method = methods{m};
        results(c).(method).TopCIvol = TopCIvol_all.(method);
        results(c).(method).TopCIsz  = TopCIsz_all.(method);
        results(c).(method).TopTF    = TopTF_all.(method);
        results(c).(method).PaS      = PaS_all.(method);
        results(c).(method).Conductance = Phi_all.(method);
        results(c).(method).NMI      = NMI_all.(method);
        results(c).(method).Fscore   = Fscore_all.(method);
        results(c).(method).time     = time_all.(method);
    end
    results(c).graph_name = graph_cases{c};
    results(c).k_est = k_est_all;
end

%% Print results table
fprintf('%-30s | %-6s | %-10s | %-9s | %-7s | %-12s | %-8s | %-9s | %-5s | %-7s\n', ...
    'Graph', 'k', 'TopCIvol', 'TopCIsz', 'TopTF', 'Conductance', 'NMI', 'F-score', 'PaS', 'Time');
fprintf(repmat('-', 1, 120)); fprintf('\n');
for m = 1:length(methods)
    method = methods{m};
   fprintf('++++++++ %s ++++++++ \n',upper(method));
    for c = 1:nc
        fprintf('%-30s | %6.2f | %10.4f | %9.2f | %7.1f | %12.4f | %8.4f | %9.4f | %5.4f | %7.2fs\n', ...
            [results(c).graph_name(1:end-4)], ...
            mean(results(c).k_est), ...
            mean(results(c).(method).TopCIvol), ...
            mean(results(c).(method).TopCIsz), ...
            mean(results(c).(method).TopTF), ...
            mean(results(c).(method).Conductance), ...
            mean(results(c).(method).NMI), ...
            mean(results(c).(method).Fscore), ...
            mean(results(c).(method).PaS), ...
            mean(results(c).(method).time));
    end
end
