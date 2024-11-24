function [x_inferred, inferred_label] = label_data(x,labels,method)
% label_data - This function reassigns predicted cluster labels to align 
%              with ground-truth labels using two approaches:
%               (1) Mode-based
%               (2) Cost-based
%
%% Syntax:
%        label_data(x, labels, method)
%
%% Input Arguments:
%       *Required Input Arguments*
%       - x:                Predicted cluster labels for the data points.
%       - labels:           True labels (ground truth) for the data points.
%       - method:           Specifies which method to use for inferring 
%                           labels: 1 for mode-based, others for
%                           cost-based.
%
%% Output Arguments:
%       - x_inferred:       Reassigned labels for each data point, 
%                           matching the ground-truth label distribution.
%       - inferred_label:   Mapping from predicted cluster labels (x) to 
%                           inferred ground-truth labels (labels).
%
% ------------------------------------------
%   Author: Dimosthenis Pasadakis
%   Source: https://github.com/DmsPas/Multiway-p-spectral-clustering
%

% Number of clusters
K = max(x);

L = size(unique(labels),1);

inferred_label = -1*ones(K,1);

x_inferred = x;


if method == 1
    
    for k=1:K
        index             = find(x==k);
        cur_labels        = labels(index);
        actual_label      = mode(cur_labels);
        inferred_label(k) = actual_label;
        x_inferred(index) = actual_label;
    end
    
%         if size(unique(inferred_label),1) ~= size(inferred_label,1)
%             fprintf('Clusters are merging into 1 label\n');
%             fprintf('Probably unbalanced dataset\n');
%         end
    
    
else
    
    x_labeled     = x;       
    confusion_mat = confusionmat(labels, x_labeled);
    Cost_mat      = zeros(L,K);
    
    for i = 1:L
        for j = 1:K
            %             delta = confusion_mat(i,:)-confusion_mat(i,j);
            Cost_mat(i,j) = sum(confusion_mat(i,:))-confusion_mat(i,j);
        end
    end
    costofnonassignment = 2*max(max(Cost_mat));
%     costofnonassignment = .2;
    [assignments, unassignedTracks, unassignedDetections] = ...
        assignmunkres(Cost_mat,costofnonassignment);
    
    
    %% Munkres vectorized implementation
    %     [assignments_Matrix,Cost_assgn]    = munkres(Cost_mat);
    %     assignments = zeros(K,2);
    %     for i = 1:K
    %         assignments(i,1) = i;
    %         assignments(i,2) = find(assignments_Matrix(i,:));
    %     end
    %%
    
    inferred_label(assignments(:,2)) = assignments(:,1);
       
    for j = 1:L
        index             = find(x_labeled==j);
        x_inferred(index) = inferred_label(j);
    end
    
    if size(unique(inferred_label),1) ~= size(inferred_label,1)
        fprintf('Clusters are merging into 1 label\n');
        fprintf('Probably unbalanced dataset\n');
    end
    
    
end



end