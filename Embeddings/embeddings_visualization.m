clear; close all;
rng(1991);

% Uncomment out the graphs cases you want to analyze
graph_cases = {
    "LFR_mu10.mat",
    %"DSBM_8blocks_5000nodes_0.000000noise.mat",
    %"DSBM_5blocks_500nodes.mat"
    %"8blocks_500nodes.mat",
    %"8blocks_500nodes_unsorted.mat",
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
            [embeddings, clusters] = SVD_unscaled_embeddings(W,k);
        case "SVD_scaled"
            [embeddings, clusters] = SVD_scaled_embeddings(W,k);
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

    % Dimensionality reduction of the embeddings using tSNE
    Y = tsne(embeddings);
   
    % Plot
    figure;
    gscatter(Y(:,1),Y(:,2),labels);
    title(sprintf("Ground truth. Embeddings: %s - %s", method{1}, graph_cases{i}), "Interpreter","none");

    figure;
    gscatter(Y(:,1),Y(:,2),clusters);
    title(sprintf("Inferred labels. Embeddings: %s - %s", method{1}, graph_cases{i}), "Interpreter","none");

end