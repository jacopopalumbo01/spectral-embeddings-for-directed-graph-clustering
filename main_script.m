clear all; clc;

% Synthetically generate a block-cycle
num_blocks = 10; % number of blocks
num_nodes  = 100; % numbert of nodes
conn_prob  = 0.8; % connection probability between consecutive blocks

[W,nodes] = GenBlockCycle(num_blocks, num_nodes, conn_prob);

G = digraph(W);
% Make this into a plotting function for the 
% eigenvalues
node_color = zeros(num_nodes, 3);
for i = 1:num_nodes
    if nodes(i) == 1
        node_color(i,:) = [1 0 0];
    elseif nodes(i) == 2
        node_color(i,:) = [0 1 0];
    else
        node_color(i,:) = [0 0 1];
    end
end
h = plot(G, "NodeColor",node_color);
h.XData = nodes + rand(size(nodes,1),1);
h.YData = nodes;
for i = 1:num_nodes
    if nodes(i) == 3
        h.YData(i) = 1;
    end
end
h.NodeLabel = {};
% till here


% Compute the transition matrix
P = TransitionMatrix(W);

% Compute the eigenvalues of the transition matrix
[V, D] = eig(P);

% Plot the cycle eigenvalues
figure;
scatter(real(D), imag(D), "rx");

% Get clusters
clusters = BCS(W, P, num_blocks);