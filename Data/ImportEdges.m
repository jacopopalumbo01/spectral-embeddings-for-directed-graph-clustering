function [W] = ImportEdges(filename)
% ImportEdges - Returns the adjacency matrix associated with the provided
% edges
%
% Input:
%   - filename:         Where the file listing the edges is located
%
% Output:
%   - W:            Adjacency matrix

% Open the matrix
edges = readmatrix(filename, "FileType", "text");

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

end