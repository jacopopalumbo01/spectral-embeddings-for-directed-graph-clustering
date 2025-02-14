function [A, tau] = StochasticBlockmodel(n, k, rho, P)
% StochasticBlockmodel - Generates a Stochastic Block Model
%
% Input:
%   - n:            Number of nodes
%   - k:            Number of blocks
%   - rho:          Block membership distribution
%   - P:            Block connection probability
%
% Output:
%   - A:            Adjacency matrix (nxn)
%   - tau:          Block membership function


    % Compute membership function
    tau = randsample(1:k, n, true, rho);
    % 
    % tau = sort(tau);
    
    % Initialize A
    A = zeros(n,n);
    
    % Populate the adjacency matric
    for i = 1:n
        for j = 1:n
            if rand() <= P(tau(i),tau(j))
                A(i,j) = 1;
            end
        end
    end
    
end