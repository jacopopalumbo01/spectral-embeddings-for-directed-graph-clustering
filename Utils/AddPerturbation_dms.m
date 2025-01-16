function [W] = AddPerturbation_dms(W, pert_prob)
    % AddPerturbation - Add perturbation to a provided weighted graph
    %
    % Input:
    %   - W: Adjacency matrix (numeric, weighted or unweighted)
    %   - pert_prob: Scalar probability of adding a perturbation (new edge or weight)
    %
    % Output:
    %   - W: Updated adjacency matrix with perturbations (numeric)
    
        % Generate a random matrix the same size as W
        random_matrix = rand(size(W));
    
        % Determine where to perturb (random values <= pert_prob)
        perturbation_mask = random_matrix <= pert_prob;
    
        % Generate new random weights for all perturbations
        new_weights = rand(size(W));
    
        % Apply perturbations
        W(perturbation_mask) = new_weights(perturbation_mask);
    
        % Prevent self-loops (set diagonal entries to 0)
        W(1:size(W, 1) + 1:end) = 0;
    
end
    