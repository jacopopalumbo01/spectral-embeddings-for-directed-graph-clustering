clear; close all;warning off;
rng(1991);

%% Graphs
dataset_path = "synthetic/generated";

graph_cases = {};

graph_cases(1).name = "DSBM_2blocks";
graph_cases(1).graphs = {
    "DSBM_2blocks_5000nodes_0.000000noise_0seed.mat",
    "DSBM_2blocks_5000nodes_0.000000noise_1234seed.mat",
    "DSBM_2blocks_5000nodes_0.000000noise_123seed.mat",
    "DSBM_2blocks_5000nodes_0.000000noise_1991seed.mat",
    "DSBM_2blocks_5000nodes_0.000000noise_1seed.mat",
    "DSBM_2blocks_5000nodes_0.000000noise_2001seed.mat",
    "DSBM_2blocks_5000nodes_0.000000noise_2025seed.mat",
    "DSBM_2blocks_5000nodes_0.000000noise_2seed.mat",
    "DSBM_2blocks_5000nodes_0.000000noise_42seed.mat",
    "DSBM_2blocks_5000nodes_0.000000noise_5seed.mat",
};

graph_cases(2).name = "DSBM_3blocks";
graph_cases(2).graphs = {
    "DSBM_3blocks_5000nodes_0.000000noise_0seed.mat",
    "DSBM_3blocks_5000nodes_0.000000noise_1234seed.mat",
    "DSBM_3blocks_5000nodes_0.000000noise_123seed.mat",
    "DSBM_3blocks_5000nodes_0.000000noise_1991seed.mat",
    "DSBM_3blocks_5000nodes_0.000000noise_1seed.mat",
    "DSBM_3blocks_5000nodes_0.000000noise_2001seed.mat",
    "DSBM_3blocks_5000nodes_0.000000noise_2025seed.mat",
    "DSBM_3blocks_5000nodes_0.000000noise_2seed.mat",
    "DSBM_3blocks_5000nodes_0.000000noise_42seed.mat",
    "DSBM_3blocks_5000nodes_0.000000noise_5seed.mat",
};

graph_cases(3).name = "DSBM_4blocks";
graph_cases(3).graphs = {
    "DSBM_4blocks_5000nodes_0.000000noise_0seed.mat",
    "DSBM_4blocks_5000nodes_0.000000noise_1234seed.mat",
    "DSBM_4blocks_5000nodes_0.000000noise_123seed.mat",
    "DSBM_4blocks_5000nodes_0.000000noise_1991seed.mat",
    "DSBM_4blocks_5000nodes_0.000000noise_1seed.mat",
    "DSBM_4blocks_5000nodes_0.000000noise_2001seed.mat",
    "DSBM_4blocks_5000nodes_0.000000noise_2025seed.mat",
    "DSBM_4blocks_5000nodes_0.000000noise_2seed.mat",
    "DSBM_4blocks_5000nodes_0.000000noise_42seed.mat",
    "DSBM_4blocks_5000nodes_0.000000noise_5seed.mat",
};

graph_cases(4).name = "DSBM_5blocks";
graph_cases(4).graphs = {
    "DSBM_5blocks_5000nodes_0.000000noise_0seed.mat",
    "DSBM_5blocks_5000nodes_0.000000noise_1234seed.mat",
    "DSBM_5blocks_5000nodes_0.000000noise_123seed.mat",
    "DSBM_5blocks_5000nodes_0.000000noise_1991seed.mat",
    "DSBM_5blocks_5000nodes_0.000000noise_1seed.mat",
    "DSBM_5blocks_5000nodes_0.000000noise_2001seed.mat",
    "DSBM_5blocks_5000nodes_0.000000noise_2025seed.mat",
    "DSBM_5blocks_5000nodes_0.000000noise_2seed.mat",
    "DSBM_5blocks_5000nodes_0.000000noise_42seed.mat",
    "DSBM_5blocks_5000nodes_0.000000noise_5seed.mat",
};

graph_cases(5).name = "DSBM_6blocks";
graph_cases(5).graphs = {
    "DSBM_6blocks_5000nodes_0.000000noise_0seed.mat",
    "DSBM_6blocks_5000nodes_0.000000noise_1234seed.mat",
    "DSBM_6blocks_5000nodes_0.000000noise_123seed.mat",
    "DSBM_6blocks_5000nodes_0.000000noise_1991seed.mat",
    "DSBM_6blocks_5000nodes_0.000000noise_1seed.mat",
    "DSBM_6blocks_5000nodes_0.000000noise_2001seed.mat",
    "DSBM_6blocks_5000nodes_0.000000noise_2025seed.mat",
    "DSBM_6blocks_5000nodes_0.000000noise_2seed.mat",
    "DSBM_6blocks_5000nodes_0.000000noise_42seed.mat",
    "DSBM_6blocks_5000nodes_0.000000noise_5seed.mat",
};

graph_cases(6).name = "DSBM_7blocks";
graph_cases(6).graphs = {
    "DSBM_7blocks_5000nodes_0.000000noise_0seed.mat",
    "DSBM_7blocks_5000nodes_0.000000noise_1234seed.mat",
    "DSBM_7blocks_5000nodes_0.000000noise_123seed.mat",
    "DSBM_7blocks_5000nodes_0.000000noise_1991seed.mat",
    "DSBM_7blocks_5000nodes_0.000000noise_1seed.mat",
    "DSBM_7blocks_5000nodes_0.000000noise_2001seed.mat",
    "DSBM_7blocks_5000nodes_0.000000noise_2025seed.mat",
    "DSBM_7blocks_5000nodes_0.000000noise_2seed.mat",
    "DSBM_7blocks_5000nodes_0.000000noise_42seed.mat",
    "DSBM_7blocks_5000nodes_0.000000noise_5seed.mat",
};

graph_cases(7).name = "DSBM_8blocks";
graph_cases(7).graphs = {
    "DSBM_8blocks_5000nodes_0.000000noise_0seed.mat",
    "DSBM_8blocks_5000nodes_0.000000noise_1234seed.mat",
    "DSBM_8blocks_5000nodes_0.000000noise_123seed.mat",
    "DSBM_8blocks_5000nodes_0.000000noise_1991seed.mat",
    "DSBM_8blocks_5000nodes_0.000000noise_1seed.mat",
    "DSBM_8blocks_5000nodes_0.000000noise_2001seed.mat",
    "DSBM_8blocks_5000nodes_0.000000noise_2025seed.mat",
    "DSBM_8blocks_5000nodes_0.000000noise_2seed.mat",
    "DSBM_8blocks_5000nodes_0.000000noise_42seed.mat",
    "DSBM_8blocks_5000nodes_0.000000noise_5seed.mat",
};

graph_cases(8).name = "DSBM_9blocks";
graph_cases(8).graphs = {
    "DSBM_9blocks_5000nodes_0.000000noise_0seed.mat",
    "DSBM_9blocks_5000nodes_0.000000noise_1234seed.mat",
    "DSBM_9blocks_5000nodes_0.000000noise_123seed.mat",
    "DSBM_9blocks_5000nodes_0.000000noise_1991seed.mat",
    "DSBM_9blocks_5000nodes_0.000000noise_1seed.mat",
    "DSBM_9blocks_5000nodes_0.000000noise_2001seed.mat",
    "DSBM_9blocks_5000nodes_0.000000noise_2025seed.mat",
    "DSBM_9blocks_5000nodes_0.000000noise_2seed.mat",
    "DSBM_9blocks_5000nodes_0.000000noise_42seed.mat",
    "DSBM_9blocks_5000nodes_0.000000noise_5seed.mat",
};

graph_cases(9).name = "DSBM_10blocks";
graph_cases(9).graphs = {
    "DSBM_10blocks_5000nodes_0.000000noise_0seed.mat",
    "DSBM_10blocks_5000nodes_0.000000noise_1234seed.mat",
    "DSBM_10blocks_5000nodes_0.000000noise_123seed.mat",
    "DSBM_10blocks_5000nodes_0.000000noise_1991seed.mat",
    "DSBM_10blocks_5000nodes_0.000000noise_1seed.mat",
    "DSBM_10blocks_5000nodes_0.000000noise_2001seed.mat",
    "DSBM_10blocks_5000nodes_0.000000noise_2025seed.mat",
    "DSBM_10blocks_5000nodes_0.000000noise_2seed.mat",
    "DSBM_10blocks_5000nodes_0.000000noise_42seed.mat",
    "DSBM_10blocks_5000nodes_0.000000noise_5seed.mat",
};

%% Methods
methods = {
    "SVD_unscaled",
    "SVD_unscaled_tSNE",
    "SVD_scaled",
    "SVD_scaled_tSNE",
    "SKEW",
    "SKEW_tSNE",
    "Herm",
    "Herm_tSNE",
    "BAS",
    "BAS_tSNE",
    "BCS",
    "BCS_tSNE"
};



%% Initialize results
nc = length(graph_cases);
nm = length(methods);
results = struct;

%% Compute results
for i = 1:nc
    fprintf("------------------------------------\n");
    fprintf("|  Evaluating graphs %s  |\n", graph_cases(i).name);
    fprintf("------------------------------------\n");
    
    ng = length(graph_cases(i).graphs);

    for j = 1:ng
        graph_name = string(graph_cases(i).graphs(j));
        fprintf("Graph: %s\n", graph_name);

        % Load graph
        fprintf("   Loading graph...\n");

        loaded = load(sprintf("%s/%s", dataset_path, graph_name));
        W      = loaded.W;
        labels = loaded.labels;
        
        % Get number of labels
        k = size(unique(labels), 1);

        for m = 1:nm
            method = string(methods(m));
            fprintf("   Computing labels using %s...\n", method);

            % Get labels
            t_start = tic;
            switch methods{m}
                case "SVD_unscaled"
                    [embeddings, clusters] = SVD_unscaled(W,k);
                case "SVD_unscaled_tSNE"
                    [embeddings, clusters] = SVD_unscaled_tSNE(W,k);
                case "SVD_scaled"
                    [embeddings, clusters] = SVD_scaled(W,k);
                case "SVD_scaled_tSNE"
                    [embeddings, clusters] = SVD_scaled_tSNE(W,k);
                case "SKEW"
                    [embeddings, clusters] = SkewSymmetricClustering(W,k);
                case "SKEW_tSNE"
                    [embeddings, clusters] = SkewSymmetricClustering_tSNE(W,k);
                case "Herm"
                    [embeddings, clusters] = HermitianClustering(W,k);
                case "Herm_tSNE"
                    [embeddings, clusters] = HermitianClustering_tSNE(W,k);
                case "BAS"
                    [embeddings, clusters] = BAS(W,k);
                case "BAS_tSNE"
                    [embeddings, clusters] = BAS_tSNE(W,k);
                case "BCS"
                    [embeddings, clusters] = BCS(W,k);
                case "BCS_tSNE"
                    [embeddings, clusters] = BCS_tSNE(W,k);
                otherwise
                    fprintf("Method provided is not implemented.\n");
                    break;
            end
            time = toc(t_start);

            fprintf("   Computing metrics...\n");

            c_thres = 2*(k-1);

            Phi = Compute_Conductance(W, clusters);
            PaS = PrecedenceAlignmentScore(W, clusters);
            [CI_vals, CIsz_vals, CIvol_vals, TopCIvol, TopCIsz, TopTF] = ComputeCI(W, clusters, c_thres);
            [NMI, Fscore] = Compute_ext_metrics(labels, clusters); 
            ARI = clustereval(clusters, labels, "ari");
        
            
            % Save results
            results.(graph_cases(i).name).(method).TopCIvol(j) = TopCIvol;
            results.(graph_cases(i).name).(method).TopCIsz(j) = TopCIsz;
            results.(graph_cases(i).name).(method).TopTF(j) = TopTF;
            results.(graph_cases(i).name).(method).PaS(j) = PaS;
            results.(graph_cases(i).name).(method).Conductance(j) = Phi;
            results.(graph_cases(i).name).(method).NMI(j) = NMI;
            results.(graph_cases(i).name).(method).Fscore(j) = Fscore;
            results.(graph_cases(i).name).(method).ARI(j) = ARI;
            results.(graph_cases(i).name).(method).time(j) = time;
        end
    end
end

fprintf("---------------------------------\n");
fprintf("| Work ended. Saving results... |\n");
fprintf("---------------------------------\n");

save("experiments/results/DSBM_increasing_clusters.mat", "results");

