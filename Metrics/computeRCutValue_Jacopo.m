function RCut = computeRCutValue_Jacopo(clusters,W)
% computeRCutValue - This function computes the Ratio Cut value 
%
%% Syntax:
%        computeRCutValue(clusters, W)
%
%% Input Arguments:
%       *Required Input Arguments*
%       - clusters:         The labels of the clusters associated with each
%                           node
%       - W:                The graph adjacency matrix
%
%% Output Arguments:
%       - RCut:             The value of Ratio Cut
%    

% Number of clusters
K = max(clusters);

% Initialized RCut value
RCut = 0;

for k = 1:K
    W2  = W(clusters==k,clusters~=k);
    
    % Get number of nodes for the current partition
    num_nodes = size(W2, 1);

    if num_nodes == 0
        RCut = RCut + 0.0;
    else
        RCut = RCut + (full(sum(sum(W2))) / num_nodes);
    end
    
end

end