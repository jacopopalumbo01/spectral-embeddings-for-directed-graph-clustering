clear;clc;
rng(47);

fprintf("------------------\n");
fprintf("EcoFlorida Dataset\n");
fprintf("------------------\n");

% Import the dataset
path_edges    = sprintf("%s/data/eco-florida/florida-bay.txt", pwd);
path_clusters = sprintf("%s/data/eco-florida/Florida-bay-clusters.csv", pwd);
% Open the matrix
edges = readmatrix(path_edges, "FileType", "text");

% Get the number of edges
num_edges = size(edges,1);
% Get the number of nodes
num_nodes = max([edges(:,1); edges(:,2)]) + 1;

% Create adjacency matrix
W = zeros(num_nodes);

% Populate the adjacency matrix
for i = 1:num_edges
    W(edges(i,1) + 1, edges(i,2) + 1) = 1;
end

T = readtable(path_clusters,'PreserveVariableNames',true);

% Extract labels
labels = -1 * ones(num_nodes,1);
for row = 1:size(T,1)
    if ~isnan(T(row,:).cluster)
        labels(row) = T(row,:).cluster;
    end
end   

% Clean unassigned nodes
unassigned = labels == -1;
W(unassigned, :) = [];
W(:,unassigned)  = [];
labels(unassigned) = [];

% Compute final values
num_nodes = size(W,1);
num_edges = sum(W,"all");

fprintf("Num nodes: %d, Num edges: %d\n", num_nodes, num_edges);

save(sprintf("%s/data/eco-florida/eco.mat", pwd), "W", "labels");