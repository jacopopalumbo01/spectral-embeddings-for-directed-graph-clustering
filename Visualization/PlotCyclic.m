function PlotCyclic(A, num_blocks, nodes)
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
%

num_nodes = size(nodes,1);

%% Generate colors for different blocks
available_colors = [255 0 0
                    0 255 0
                    0 0 255
                    255 0 255
                    0 255 255] / 255;


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
plot(G, "NodeColor", node_color);

h.XData = nodes_pos(:,1) + rand(size(nodes,1),1)/10;
h.YData = nodes_pos(:,2) + rand(size(nodes,1),1)/10;

h.NodeLabel = {};
end

