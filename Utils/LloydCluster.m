function [membership, centroids] = LloydCluster(X, k, max_iters, tol)
% LloydCluster - Lloyd’s Algorithm
%
%% Syntax:
%        LloydCluster(X, k, max_iters, tol)
%
%% Input Arguments:
%       *Required Input Arguments*
%       - X:                Matrix. Each row is a point of the set on which
%                           the clustering will be performed.
%       - k:                Number of clusters
%       - max_iters:        Maximum number of iterations
%       - tol:              Minimum tolerance
%
%% Output Arguments:
%       - membership:       Vector with mappings (point -> cluster)
%       - centroids:        The centroids found by the algorithm
%

    if nargin == 2
        max_iters = 100;
        tol = 1e-6;
    end
    if nargin == 3
        tol = 1e-6;
    end
    

    % Get dimension
    num_points = size(X,1);
    dim = size(X, 2);
    
    %% Initialize random centroids taken from X
    centroids = X(randperm(num_points, k), :);
    new_centroids = centroids;

    % Initialize other variables
    distance = zeros(k,1);
    membership = zeros(num_points, 1);

    for iter = 1:max_iters
        %% Compute nearest centroid for each point
        for j = 1:num_points
            % Compute distance
            for c = 1:k
                sub = transpose(X(j,:)) - centroids(c);
                distance(c) = norm(sub);
            end

            [~, membership(j)] = min(distance);
        end
       

        %% Update centroids
        residual = zeros(k, 1);

        for j = 1:k
            cluster = X(membership == j,:);
            
            if ~isempty(cluster)
                new_centroids(j,:) = mean(cluster);

                % Compute residual (distance from previous step centroid
                sub = new_centroids(j) - centroids(j);
                residual(j) = norm(sub);
            else
                % Reinitialize empty cluster centroid to a random point
                new_centroids(j,:) = X(randi(num_points),:);
            end

        end

        % Update
        centroids = new_centroids;
        
        % Check stop condition
        if max(residual) <= tol
            fprintf('Lloyd''s Algorithm converged at iteration %d\n', iter);
            break;
        end
    end

    if iter == max_iters
        fprintf('Lloyd''s Algorithm did not reach the required tol (%d) at iter=%d\n', tol, max_iters);
    end

end