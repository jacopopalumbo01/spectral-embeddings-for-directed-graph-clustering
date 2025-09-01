clear; close all;
rng(1991);

% Uncomment out the graphs cases you want to analyze
graph_cases = {
    "LFR-2-0.10.mat",
    %"DSBM_3blocks_5000nodes_0.000000noise_0seed.mat",
    %"DSBM_4blocks_5000nodes_0.000000noise_0seed.mat",
    %"DSBM_5blocks_5000nodes_0.000000noise_0seed.mat",
    %"DSBM_10blocks_5000nodes_0.000000noise_0seed.mat",
};

% Uncomment out ONLY ONE method you want to use
method = {
  % "SVD_unscaled",
   "SVD_scaled",
  % "BAS",
  % "BAS_tSNE",
  % "BCS",
  % "BCS_tSNE"
};

%% Initialize results
nc = length(graph_cases);

fprintf("-------------------------------\n");
fprintf("Method: %s\n", method{1});
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
    switch method{1}
        case "SVD_unscaled"
            [embeddings, clusters] = SVD_unscaled(W,k);
        case "SVD_scaled"
            [embeddings, clusters] = SVD_scaled(W,k);
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

    % Align inferred labels
    [clusters,~] = label_data(clusters,labels,2);


    % Dimensionality reduction of the embeddings using tSNE
    Y = tsne(embeddings);
   
    % Plot
    figure;
    gscatter(Y(:,1),Y(:,2),labels,[],[],25);
    title("Ground truth", "Interpreter","latex");
    legend("boxoff");
    xlabel("$x$", "Interpreter","latex");
    ylabel("$y$", "Interpreter","latex");
    set(gca, "fontsize", 50);
    if graph_cases{i} == "LFR-2-0.10.mat"
        legend("hide");
    end

    figure;
    gscatter(Y(:,1),Y(:,2),clusters,[],[],25);
    title("Inferred labels", "Interpreter","latex");
    legend("boxoff");
    xlabel("$x$", "Interpreter","latex");
    ylabel("$y$", "Interpreter","latex");
    set(gca, "fontsize", 50);
    if graph_cases{i} == "LFR-2-0.10.mat"
        legend("hide");
    end

    % Compute accuracy
    [NMI, Fscore] = Compute_ext_metrics(labels, clusters); 
    fprintf("Accuracy (F-Score): %f\n", Fscore);
end