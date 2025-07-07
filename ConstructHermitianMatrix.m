function He = ConstructHermitianMatrix(W)
% Construct the Hermitian matrix H as in https://ieeexplore.ieee.org/abstract/document/10020413
% W is the adjacency matrix (n x n), assumed binary or weighted
% Output H is a complex Hermitian matrix

    n = size(W, 1);
    if ~issparse(W) 
        He = complex(zeros(n, n));
    else
        He = sparse(n);
    end
    
    for u = 1:n
        for v = 1:n
            if W(u,v) > 0 && W(v,u) == 0
                He(u,v) = 1i * W(u,v);
            elseif W(u,v) == 0 && W(v,u) > 0
                He(u,v) = -1i * W(v,u);
            elseif W(u,v) > 0 && W(v,u) > 0
                % Replace reciprocal edges with directed difference
                if W(u,v) >= W(v,u)
                    He(u,v) = 1i * (W(u,v) - W(v,u));
                else
                    He(u,v) = -1i * (W(v,u) - W(u,v));
                end
            end
        end
    end

end
