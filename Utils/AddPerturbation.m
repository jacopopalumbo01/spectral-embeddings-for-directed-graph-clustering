function [W] = AddPerturbation(W, pert_prob)
% AddPerturbation - Add perturbation to a provided unweighted graph
%% Syntax:
%        AddPerturbation(W, pert_prob)
%
%% Input Arguments:
%       *Required Input Arguments*
%       - W:                Adjacency matrix
%       - pert_prob:        Perturbation probability. It's just a scalar
%                           
%% Output:
%       - W:                Updated adjacency matrix
%

for i = 1:size(W,1)
    for j = 1:size(W,2)
        if rand(1) <= pert_prob
            W(i,j) = 1;
        end
    end
end
end