function [Q, H] = Arnoldi(A, n)
% Arnoldi - MATLAB implementation of the Arnoldi iteration 
%
%% Syntax:
%        [Q, H]     = Arnoldi(A, n)
%
%% Description:
%       See <a href="https://en.wikipedia.org/wiki/Arnoldi_iteration">Arnoldi</a>
%
%% Input Arguments:
%       *Required Input Arguments*
%       - A:    Input Matrix (m x m)
%       - n:    im(K) - 1 = degree(K), where K is the Krylov subspace
%
%% Output Arguments:
%       - Q:    Matrix (m x (n + 1)) where the columns are orthonormal basis of the Krylov
%       - H:    Matrix ((n + 1) x n). A on basis Q. It is an upper Hessenberg matrix.
%
% ------------------------------------------
%   Author: Jacopo Palumbo
%   Email:  jacopo.palumbo@usi.ch

m = size(A, 1);

% Generate the initial arbitrary vector
q = ones(m,1);
% Normalize the vector
q = q ./ norm(q);

% Initialize matrix Q
Q = zeros(m, n + 1);
Q(:,1) = q;

% Initialize matrix h
H = zeros(n + 1, n);

for k = 2:n
    % Compute the new candidate vector
    Q(:,k) = A * Q(:, k - 1);
    
    % Substract the projections on previous vectors
    for j = 1:k - 1
        H(j, k - 1) = dot(conj(Q(:, j)), Q(:, k));
        Q(:, k) = Q(:, k) - Q(:, j) .* H(j, k - 1);
    end
    H(k, k - 1) = norm(Q(:, k));
    % Add the produced vector to Q
    Q(:, k) = Q(:, k) ./ H(k, k - 1);
end

end