function [W, labels] = GenerateGraph(n,k,rho,P,sorted,epsilon)
% GenerateGraph - Generates a Stochastic Block Model graph
%
% Input:
%   - n:            Number of nodes
%   - k:            Number of blocks
%   - rho:          Block membership distribution
%   - P:            Block connection probability
%   - sorted:       Sort the adjacency matrix based on the nodes order
%   - epsilon:      *optional* Perturbation magnitude
%
% Output:
%   - W:            Adjacency matrix (nxn)
%   - labels:       Block membership function

if nargin < 5; sorted = 0; end
if nargin < 6; epsilon = 0; end

if epsilon > 1
    error("Epsilon should be less than 1. %f was provided instead", epsilon);
end

% Generate unperturbed graph
[W, labels] = StochasticBlockmodel(n,k,rho,P,sorted);

% Add perturbation
if epsilon ~= 0
    % Create perturbing block connection probability matrix
    Q  = rand(k,k);

    % Computing new perturbing block connection probability matrix given
    % the magnitude of perturbation epsilon
    P_pert = epsilon * Q;

    % Generate perturbing graph
    [W_pert, ~] = StochasticBlockmodel(n,k,rho,P_pert);

    % Combine unperturbed and perturbing graph
    W = CombineBlockmodels(W,W_pert);
end
    
end