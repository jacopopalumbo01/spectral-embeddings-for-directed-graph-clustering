clear; close all;warning off;
rng(1991);

%% Graphs
dataset_path = "synthetic/LFR/generated";

graph_cases = {};

graph_cases(1).name = "LFR_10";
graph_cases(1).graphs = {
    "LFR-0-0.10.mat",
    "LFR-1234-0.10.mat",
    "LFR-123-0.10.mat",
    "LFR-1991-0.10.mat",
    "LFR-1-0.10.mat",
    "LFR-2001-0.10.mat",
    "LFR-2025-0.10.mat",
    "LFR-2-0.10.mat",
    "LFR-42-0.10.mat",
    "LFR-5-0.10.mat",
};

graph_cases(2).name = "LFR_12";
graph_cases(2).graphs = {
    "LFR-0-0.12.mat",
    "LFR-1234-0.12.mat",
    "LFR-123-0.12.mat",
    "LFR-1991-0.12.mat",
    "LFR-1-0.12.mat",
    "LFR-2001-0.12.mat",
    "LFR-2025-0.12.mat",
    "LFR-2-0.12.mat",
    "LFR-42-0.12.mat",
    "LFR-5-0.12.mat",
};

graph_cases(3).name = "LFR_14";
graph_cases(3).graphs = {
    "LFR-0-0.14.mat",
    "LFR-1234-0.14.mat",
    "LFR-123-0.14.mat",
    "LFR-1991-0.14.mat",
    "LFR-1-0.14.mat",
    "LFR-2001-0.14.mat",
    "LFR-2025-0.14.mat",
    "LFR-2-0.14.mat",
    "LFR-42-0.14.mat",
    "LFR-5-0.14.mat",
};

graph_cases(4).name = "LFR_16";
graph_cases(4).graphs = {
    "LFR-0-0.16.mat",
    "LFR-1234-0.16.mat",
    "LFR-123-0.16.mat",
    "LFR-1991-0.16.mat",
    "LFR-1-0.16.mat",
    "LFR-2001-0.16.mat",
    "LFR-2025-0.16.mat",
    "LFR-2-0.16.mat",
    "LFR-42-0.16.mat",
    "LFR-5-0.16.mat",
};

graph_cases(5).name = "LFR_18";
graph_cases(5).graphs = {
    "LFR-0-0.18.mat",
    "LFR-1234-0.18.mat",
    "LFR-123-0.18.mat",
    "LFR-1991-0.18.mat",
    "LFR-1-0.18.mat",
    "LFR-2001-0.18.mat",
    "LFR-2025-0.18.mat",
    "LFR-2-0.18.mat",
    "LFR-42-0.18.mat",
    "LFR-5-0.18.mat",
};

graph_cases(6).name = "LFR_20";
graph_cases(6).graphs = {
    "LFR-0-0.20.mat",
    "LFR-1234-0.20.mat",
    "LFR-123-0.20.mat",
    "LFR-1991-0.20.mat",
    "LFR-1-0.20.mat",
    "LFR-2001-0.20.mat",
    "LFR-2025-0.20.mat",
    "LFR-2-0.20.mat",
    "LFR-42-0.20.mat",
    "LFR-5-0.20.mat",
};

graph_cases(7).name = "LFR_22";
graph_cases(7).graphs = {
    "LFR-0-0.22.mat",
    "LFR-1234-0.22.mat",
    "LFR-123-0.22.mat",
    "LFR-1991-0.22.mat",
    "LFR-1-0.22.mat",
    "LFR-2001-0.22.mat",
    "LFR-2025-0.22.mat",
    "LFR-2-0.22.mat",
    "LFR-42-0.22.mat",
    "LFR-5-0.22.mat",
};

graph_cases(8).name = "LFR_24";
graph_cases(8).graphs = {
    "LFR-0-0.24.mat",
    "LFR-1234-0.24.mat",
    "LFR-123-0.24.mat",
    "LFR-1991-0.24.mat",
    "LFR-1-0.24.mat",
    "LFR-2001-0.24.mat",
    "LFR-2025-0.24.mat",
    "LFR-2-0.24.mat",
    "LFR-42-0.24.mat",
    "LFR-5-0.24.mat",
};

graph_cases(9).name = "LFR_26";
graph_cases(9).graphs = {
    "LFR-0-0.26.mat",
    "LFR-1234-0.26.mat",
    "LFR-123-0.26.mat",
    "LFR-1991-0.26.mat",
    "LFR-1-0.26.mat",
    "LFR-2001-0.26.mat",
    "LFR-2025-0.26.mat",
    "LFR-2-0.26.mat",
    "LFR-42-0.26.mat",
    "LFR-5-0.26.mat",
};

graph_cases(10).name = "LFR_28";
graph_cases(10).graphs = {
    "LFR-0-0.28.mat",
    "LFR-1234-0.28.mat",
    "LFR-123-0.28.mat",
    "LFR-1991-0.28.mat",
    "LFR-1-0.28.mat",
    "LFR-2001-0.28.mat",
    "LFR-2025-0.28.mat",
    "LFR-2-0.28.mat",
    "LFR-42-0.28.mat",
    "LFR-5-0.28.mat",
};

graph_cases(11).name = "LFR_30";
graph_cases(11).graphs = {
    "LFR-0-0.30.mat",
    "LFR-1234-0.30.mat",
    "LFR-123-0.30.mat",
    "LFR-1991-0.30.mat",
    "LFR-1-0.30.mat",
    "LFR-2001-0.30.mat",
    "LFR-2025-0.30.mat",
    "LFR-2-0.30.mat",
    "LFR-42-0.30.mat",
    "LFR-5-0.30.mat",
};

graph_cases(12).name = "LFR_32";
graph_cases(12).graphs = {
    "LFR-0-0.32.mat",
    "LFR-1234-0.32.mat",
    "LFR-123-0.32.mat",
    "LFR-1991-0.32.mat",
    "LFR-1-0.32.mat",
    "LFR-2001-0.32.mat",
    "LFR-2025-0.32.mat",
    "LFR-2-0.32.mat",
    "LFR-42-0.32.mat",
    "LFR-5-0.32.mat",
};

graph_cases(13).name = "LFR_34";
graph_cases(13).graphs = {
    "LFR-0-0.34.mat",
    "LFR-1234-0.34.mat",
    "LFR-123-0.34.mat",
    "LFR-1991-0.34.mat",
    "LFR-1-0.34.mat",
    "LFR-2001-0.34.mat",
    "LFR-2025-0.34.mat",
    "LFR-2-0.34.mat",
    "LFR-42-0.34.mat",
    "LFR-5-0.34.mat",
};

graph_cases(14).name = "LFR_36";
graph_cases(14).graphs = {
    "LFR-0-0.36.mat",
    "LFR-1234-0.36.mat",
    "LFR-123-0.36.mat",
    "LFR-1991-0.36.mat",
    "LFR-1-0.36.mat",
    "LFR-2001-0.36.mat",
    "LFR-2025-0.36.mat",
    "LFR-2-0.36.mat",
    "LFR-42-0.36.mat",
    "LFR-5-0.36.mat",
};

graph_cases(15).name = "LFR_38";
graph_cases(15).graphs = {
    "LFR-0-0.38.mat",
    "LFR-1234-0.38.mat",
    "LFR-123-0.38.mat",
    "LFR-1991-0.38.mat",
    "LFR-1-0.38.mat",
    "LFR-2001-0.38.mat",
    "LFR-2025-0.38.mat",
    "LFR-2-0.38.mat",
    "LFR-42-0.38.mat",
    "LFR-5-0.38.mat",
};

graph_cases(16).name = "LFR_40";
graph_cases(16).graphs = {
    "LFR-0-0.40.mat",
    "LFR-1234-0.40.mat",
    "LFR-123-0.40.mat",
    "LFR-1991-0.40.mat",
    "LFR-1-0.40.mat",
    "LFR-2001-0.40.mat",
    "LFR-2025-0.40.mat",
    "LFR-2-0.40.mat",
    "LFR-42-0.40.mat",
    "LFR-5-0.40.mat",
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

save("experiments/results/LFR_increasing_noise.mat", "results");

