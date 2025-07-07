% Parameters
n = 1000;             % number of points
k = 4;                % number of clusters
noise = 0.05;         % noise level
markerSize = 50;      % size of scatter points

% Generate Swiss roll
t = (3 * pi / 2) * (1 + 2 * rand(1, n));   % unrolled coordinate
height = 21 * rand(1, n);
x = t .* cos(t);
y = height;
z = t .* sin(t);

% Add noise
x = x + noise * randn(1, n);
y = y + noise * randn(1, n);
z = z + noise * randn(1, n);

% Assign clusters based on t
edges = linspace(min(t), max(t), k+1);
cluster_labels = discretize(t, edges);

% Stack data into matrix
X = [x; y; z]';

% --- Plot 3D Swiss Roll with clusters ---
figure;
scatter3(x, y, z, markerSize, cluster_labels, 'filled');
xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('3D Swiss Roll with %d Clusters', k));
axis equal; grid on;
view(3); colorbar;

% --- Apply t-SNE to reduce to 2D ---
Y = tsne(X, 'NumDimensions', 2, 'Perplexity', 30, 'Standardize', true);

% --- Plot 2D t-SNE embedding ---
figure;
gscatter(Y(:,1), Y(:,2), cluster_labels, [], [], 25);
title(sprintf('2D t-SNE Embedding of Swiss Roll (%d Clusters)', k));
xlabel('t-SNE 1'); ylabel('t-SNE 2');
axis equal; grid on;
