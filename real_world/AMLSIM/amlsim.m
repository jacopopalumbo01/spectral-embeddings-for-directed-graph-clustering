path_edges    = "real_world/AMLSIM/transactions_100.csv";
path_clusters = "real_world/AMLSIM/accounts_100.csv";

% open tables
edge_T = readtable(path_edges,'PreserveVariableNames',true);
clusters_T = readtable(path_clusters,'PreserveVariableNames',true);

% Number of nodes
n = size(clusters_T,1);

% Get number of edges
num_edges = size(edge_T,1);


% Create adjacency matrix
W = zeros(n);

% Populate the adjacency matrix
for i = 1:num_edges
    W(edge_T(i,:).SENDER_ACCOUNT_ID + 1, edge_T(i,:).RECEIVER_ACCOUNT_ID + 1) = 1;
end

% Save
save("real_world/AMLSIM/AMLSIM-100.mat", "W");

% Next
path_edges    = "real_world/AMLSIM/transactions_1k.csv";
path_clusters = "real_world/AMLSIM/accounts_1k.csv";

% open tables
edge_T = readtable(path_edges,'PreserveVariableNames',true);
clusters_T = readtable(path_clusters,'PreserveVariableNames',true);

% Number of nodes
n = size(clusters_T,1);

% Get number of edges
num_edges = size(edge_T,1);


% Create adjacency matrix
W = zeros(n);

% Populate the adjacency matrix
for i = 1:num_edges
    W(edge_T(i,:).SENDER_ACCOUNT_ID + 1, edge_T(i,:).RECEIVER_ACCOUNT_ID + 1) = 1;
end

% Save
save("real_world/AMLSIM/AMLSIM-1k.mat", "W");