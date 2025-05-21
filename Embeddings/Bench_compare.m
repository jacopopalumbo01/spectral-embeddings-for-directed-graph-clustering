clear; close all;warning on;
rng(1991);

% Uncomment out the graphs cases you want to analyze
graph_cases = {
    "DSBM_5blocks_5000nodes.mat",
    %"8blocks_500nodes.mat",
    %"8blocks_500nodes_unsorted.mat",
    %"10blocks_1000nodes.mat",
    %"10blocks_100nodes.mat",
    %"10blocks_200nodes.mat",
    %"10blocks_500nodes.mat",
    %"2blocks_1000nodes.mat",
    %"2blocks_100nodes.mat",
    %"2blocks_200nodes.mat",
    %"2blocks_500nodes.mat",
    %"3blocks_1000nodes.mat",
    %"3blocks_100nodes.mat",
    %"3blocks_200nodes.mat",
    %"3blocks_500nodes.mat",
    %"4blocks_1000nodes.mat",
    %"4blocks_100nodes.mat",
    %"4blocks_200nodes.mat",
    %"4blocks_500nodes.mat",
    %"5blocks_1000nodes.mat",
    %"5blocks_100nodes.mat",
    %"5blocks_200nodes.mat",
    %"5blocks_500nodes.mat",
};

% Uncomment out the methods you want to use
methods = {
    "SVD_unscaled",
    "SVD_scaled",
    "SKEW",
    "BAS",
    "BAS_tSNE",
    "BCS",
    "BCS_tSNE"
};

%% Initialize results
nc = length(graph_cases);
nm = length(methods);
results = struct;

for m = 1:nm
fprintf("-------------------------------\n");
fprintf("Method: %s\n", methods{m});
fprintf("-------------------------------\n");

    for i = 1:nc
        fprintf("Evaluating graph: %s\n", graph_cases{i});
        % Load graph
        loaded = load(graph_cases{i});
        W      = loaded.W;
        labels = loaded.labels;
        
        % Get number of labels
        k = size(unique(labels), 1);
    
        % Get embeddings
        t_start = tic;
        switch methods{m}
            case "SVD_unscaled"
                [embeddings, clusters] = SVD_unscaled_embeddings(W,k);
            case "SVD_scaled"
                [embeddings, clusters] = SVD_scaled_embeddings(W,k);
            case "SKEW"
                [embeddings, clusters] = SkewSymmetricClustering_embeddings(W,k);
            case "BAS"
                [embeddings, clusters] = BAS_embeddings(W,k);
            case "BAS_tSNE"
                [embeddings, clusters] = BAS_tSNE_embeddings(W,k);
            case "BCS"
                [embeddings, clusters] = BCS_embeddings(W,k);
            case "BCS_tSNE"
                [embeddings, clusters] = BCS_tSNE_embeddings(W,k);
            otherwise
                fprintf("Method provided is not implemented.\n");
                break;
        end
        time = toc(t_start);

        % Evaluate metrics
        c_thres = 2*(k-1);

        Phi = Compute_Conductance(W, clusters);
        PaS = PrecedenceAlignmentScore(W, clusters);
        [CI_vals, CIsz_vals, CIvol_vals, TopCIvol, TopCIsz, TopTF] = ComputeCI(W, clusters, c_thres);
        [NMI, Fscore] = Compute_ext_metrics(labels, clusters); 

        % Save metrics
        method = methods{m};

        results(i).(method).TopCIvol = TopCIvol;
        results(i).(method).TopCIsz  = TopCIsz;
        results(i).(method).TopTF    = TopTF;
        results(i).(method).PaS      = PaS;
        results(i).(method).Conductance = Phi;
        results(i).(method).NMI      = NMI;
        results(i).(method).Fscore   = Fscore;
        results(i).(method).time     = time;
    end
end

%% Print results table
fprintf('%-30s | %-10s | %-9s | %-7s | %-12s | %-8s | %-9s | %-5s | %-7s\n', ...
    'Graph', 'TopCIvol', 'TopCIsz', 'TopTF', 'Conductance', 'NMI', 'F-score', 'PaS', 'Time');
fprintf(repmat('-', 1, 120)); fprintf('\n');
for m = 1:nm
    method = methods{m};
    fprintf('++++++++ %s ++++++++ \n', upper(method));
    for c = 1:nc
        fprintf('%-30s | %10.2f | %9.4f | %7.2f | %12.1f | %8.4f | %9.4f | %5.4f | %7.2f\n', ...
            graph_cases{c}, ...
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