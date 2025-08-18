clear; close all;

% Path to raw graphs and save path
raw_path  = "synthetic/LFR/raw_graphs";
save_path = "synthetic/LFR/generated";

% Seeds and mus
seeds = [0; 1; 2; 5; 42; 123; 1234; 1991; 2001; 2025];
mus   = (0.10:0.02:0.40)';

for s = 1:size(seeds,1)
    for m = 1:size(mus,1)
        seed = seeds(s);
        mu   = mus(m);
        
        fprintf("Importing LFR with seed=%d and mu=%0.2f...\n", seed, mu);
        % Import edges and labels
        edges = importdata(sprintf("%s/network-%d-%0.2f.dat", raw_path, seed, mu));
        labels = importdata(sprintf("%s/community-%d-%0.2f.dat", raw_path, seed, mu));
        labels = labels(:,2);
        
        % Get number of nodes
        num_nodes = size(labels,1);
        
        % Initialize adjacency
        W = zeros(num_nodes);
        
        for i = 1:size(edges,1)
            W(edges(i,1), edges(i,2)) = 1;
        end
        
        fprintf("    Saving...\n");
        save(sprintf("%s/LFR-%d-%0.2f.mat", save_path, seed, mu), "W", "labels");
    end
end