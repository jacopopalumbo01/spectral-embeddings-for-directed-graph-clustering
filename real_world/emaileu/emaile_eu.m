clear;close all;
rng(47);

fprintf("------------------\n");
fprintf(" Email EU Dataset\n");
fprintf("------------------\n");

% Import the dataset
path_edges    = "real_world/emaileu/email-Eu-core.txt";
path_clusters = "real_world/emaileu/email-Eu-core-department-labels.txt";

% Open the matrix
edges = readmatrix(path_edges, "FileType", "text");
% Fix indeces to start from 1
edges = edges + ones(size(edges));

% Get the number of edges
num_edges = size(edges,1);
% Get the number of nodes
num_nodes = max([edges(:,1); edges(:,2)]);

% Create adjacency matrix
W = zeros(num_nodes);

% Populate the adjacency matrix
for i = 1:num_edges
    W(edges(i,1), edges(i,2)) = 1;
end

% Extract labels
labels = readmatrix(path_clusters, "FileType","text");
labels(:,1) = [];
% Fix indeces to start from 1
labels = labels + ones(size(labels));

% Compute final values
num_nodes = size(W,1);
num_edges = sum(W,"all");

fprintf("Num nodes: %d, Num edges: %d\n", num_nodes, num_edges);

save(sprintf("%s/data/emaileu/emaileu.mat", pwd), "W", "labels");

% Extract emailEu12 and emailEu23
[clusters_size, unique_clusters] = hist(labels,unique(labels));

% Identify 3 biggest clusters (departments)
[unique_labels, ~, label_ids] = unique(labels);
cluster_sizes = accumarray(label_ids, 1);
[~, sorted_indices] = sort(cluster_sizes, 'descend');
biggest_clusters = unique_labels(sorted_indices(1:3));

% Extract subgraph of departments 1 and 2
labels_12_mask = (labels == biggest_clusters(1)) | (labels == biggest_clusters(2));
selected_idx = find(labels_12_mask);
W_12 = W(selected_idx, selected_idx);
labels_12 = labels(selected_idx);

% Extract subgraph of departments 2 and 3
labels_23_mask = (labels == biggest_clusters(2)) | (labels == biggest_clusters(3));
selected_idx = find(labels_23_mask);
W_23 = W(selected_idx, selected_idx);
labels_23 = labels(selected_idx);

% Save them
W = W_12;
labels = labels_12;
save(sprintf("%s/data/emaileu/emaileu12.mat", pwd), "W", "labels");

W = W_23;
labels = labels_23;
save(sprintf("%s/data/emaileu/emaileu23.mat", pwd), "W", "labels");