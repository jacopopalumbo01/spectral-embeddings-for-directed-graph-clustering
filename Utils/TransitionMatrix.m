function [P] = TransitionMatrix(W)
    num_nodes = size(W, 1);
    P = W;
    for i = 1:num_nodes
        out_degree = sum(P(i,:));
    
        if out_degree ~= 0
            P(i,:) = P(i,:) ./ out_degree;
        else
            P(i,:) = 0;
        end
    end
end