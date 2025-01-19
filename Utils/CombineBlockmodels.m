function [A] = CombineBlockmodels(A1,A2)
% CombineBlockmodels - Combines two blockmodels
%
% Input:
%   - A1:           First Blockmodel
%   - A2:           Second Blockmodel
%
% Output:
%   - A:            Adjacency matrix (nxn)

    % Union of the edges
    A = A1 + A2;
    
    % Since unweighted we have to normalize all weights to 1
    positive_mask = A > 0;
    A(positive_mask) = 1;

end