function PlotCyclic(A, num_blocks, nodes, graph_name)
% PlotCyclic - Plots a Block-Cyclic Graph
%
%% Syntax:
%        PlotCyclic(A, num_blocks, nodes)
%
%% Input Arguments:
%       *Required Input Arguments*
%       - A:                Adjacency Matrix
%       - num_blocks:       Number of blocks
%       - nodes:            Nodes labels
%       - graph_name:       Graph name used during visualization

% Set default graph name
if nargin < 4
    graph_name = "";
end

num_nodes = size(nodes,1);

%% Generate colors for different blocks
available_colors = [0 0.4470 0.7410
    0.8500 0.3250 0.0980
    0.9290 0.6940 0.1250
    0.4940 0.1840 0.5560
    0.4660 0.6740 0.1880];


colors = zeros(num_blocks,3);
for i = 1:num_blocks
    colors(i,:) = available_colors(size(available_colors,1) - ...
        rem(i, size(available_colors,1)),:);
end

%% Assign colors to nodes
node_color = zeros(num_nodes, 3);

for i = 1:num_nodes
    node_color(i,:) = colors(nodes(i),:);
end

%% Get position for each block
% We get the position by diving a circle of radius 1 in num_blocks parts
degree = 360 / num_blocks;
radius = 5;

block_pos = zeros(num_blocks,2);

for i = 1:num_blocks
    block_pos(i,1) = radius*cosd(degree * i);
    block_pos(i,2) = radius*sind(degree * i);
end

%% Assign position to each node
nodes_pos = zeros(size(nodes, 1),2);
for i = 1:num_nodes
    nodes_pos(i,:) = block_pos(nodes(i),:);
end
%% Plot
G = digraph(A);

figure;
h = plot(G, "NodeColor", node_color, "MarkerSize",10, "EdgeColor","black", "LineWidth",0.5);
hold on;
h.XData = nodes_pos(:,1) + rand(size(nodes,1),1)*1.5;
h.YData = nodes_pos(:,2) + rand(size(nodes,1),1)*1.5;

axis off;

title(graph_name);

h.NodeLabel = {};

%% Plot adjacency matrix sparsity pattern
figure;
axis on;
spy(A, 'k.', 15);
axis off;
if ~isempty(graph_name)
    title(sprintf("%s Adjacency Matrix", graph_name));
end

end

