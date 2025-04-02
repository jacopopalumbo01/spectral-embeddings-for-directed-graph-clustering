function [P] = ConnectionProbabilityMatrix(type, k, conn_prob)
% ConnectionProbabilityMatrix - Returns the connection probability matrix
%
% Input:
%   - type:         Either "cyclic", "nested" or "acyclic"
%   - k:            Number of blocks
%   - conn_prob:    The probability to be used in the matrix
%
% Output:
%   - P:            The connection probability matrix

P = zeros(k,k);

if type == "cyclic"
    for i = 1:k
        for j = 1:k
            if i + 1 == j || (i == k && j == 1)
                P(i,j) = conn_prob;
            end
        end
    end
elseif type == "nested"
    for i = 1:k
        for j = 1:k
            if i < j || i == k
                P(i,j) = conn_prob;
            end
        end
    end
elseif type == "acyclic"
    for i = 1:k
        for j = 1:k
            if i < j
                P(i,j) = conn_prob;
            end
        end
    end
else 
    error("Type should be 'cyclic' or 'acyclic'. Instead '%s' was provided.", type);
end

end