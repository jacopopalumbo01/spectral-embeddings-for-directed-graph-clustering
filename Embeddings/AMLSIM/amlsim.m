path_edges    = sprintf("%s/Embeddings/AMLSIM/transactions.csv", pwd);
path_clusters = sprintf("%s/Embeddings/AMLSIM/accounts.csv", pwd);

% open tables
edge_T = readtable(path_edges,'PreserveVariableNames',true);
clusters_T = readtable(path_clusters,'PreserveVariableNames',true);

% Number of nodes is fixed to 100000
n = size(clusters_T,1);

% Get number of edges
num_edges = size(edge_T,1);


% Create adjacency matrix
W = zeros(n);

% Populate the adjacency matrix
for i = 1:num_edges
    W(edge_T(i,:).SENDER_ACCOUNT_ID + 1, edge_T(i,:).RECEIVER_ACCOUNT_ID + 1) = 1;
end


% Get unique banks
%banks = unique(clusters_T.bank_id);

% Extract labels
%labels = -1 * ones(n,1);
%for row = 1:n
%    labels(row) = find(strcmp(banks, clusters_T(row,:).bank_id));
%end   