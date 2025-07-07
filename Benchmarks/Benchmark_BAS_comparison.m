clear; clc; close all;
rng(47);

addpaths_SC(); % Make sure your BAS functions are in path

%% Parameters
graph_cases = {
    %'BlockAcyclicPerturbation/3-blocks_2500-nodes_0.000000-pert';
    %'BlockAcyclicPerturbation/3-blocks_2500-nodes_0.100000-pert';
    %'BlockAcyclicPerturbation/3-blocks_2500-nodes_0.200000-pert';
    %'BlockAcyclicPerturbation/3-blocks_2500-nodes_0.300000-pert';
    %'BlockAcyclicPerturbation/3-blocks_2500-nodes_0.400000-pert';
    %'BlockAcyclicPerturbation/3-blocks_2500-nodes_0.500000-pert';
    %'BlockAcyclicPerturbation/3-blocks_2500-nodes_0.600000-pert';
    %'BlockAcyclicPerturbation/3-blocks_2500-nodes_0.700000-pert';
    %'BlockAcyclicPerturbation/3-blocks_2500-nodes_0.800000-pert'; 
    %'BlockAcyclicPerturbation/3-blocks_2500-nodes_0.900000-pert';
    %'BlockAcyclicPerturbation/5-blocks_2500-nodes_0.000000-pert';
    %'BlockAcyclicPerturbation/5-blocks_2500-nodes_0.100000-pert';
    %'BlockAcyclicPerturbation/5-blocks_2500-nodes_0.200000-pert';
    %'BlockAcyclicPerturbation/5-blocks_2500-nodes_0.300000-pert';
    %'BlockAcyclicPerturbation/5-blocks_2500-nodes_0.400000-pert';
    %'BlockAcyclicPerturbation/5-blocks_2500-nodes_0.500000-pert';
    %'BlockAcyclicPerturbation/5-blocks_2500-nodes_0.600000-pert';
    %'BlockAcyclicPerturbation/5-blocks_2500-nodes_0.700000-pert';
    %'BlockAcyclicPerturbation/5-blocks_2500-nodes_0.800000-pert';
    %'BlockAcyclicPerturbation/5-blocks_2500-nodes_0.900000-pert';
    'BlockAcyclicPerturbation/8-blocks_2500-nodes_0.000000-pert';
    'BlockAcyclicPerturbation/8-blocks_2500-nodes_0.100000-pert';
    'BlockAcyclicPerturbation/8-blocks_2500-nodes_0.200000-pert';
    'BlockAcyclicPerturbation/8-blocks_2500-nodes_0.300000-pert';
    'BlockAcyclicPerturbation/8-blocks_2500-nodes_0.400000-pert';
    'BlockAcyclicPerturbation/8-blocks_2500-nodes_0.500000-pert';
    'BlockAcyclicPerturbation/8-blocks_2500-nodes_0.600000-pert';
    'BlockAcyclicPerturbation/8-blocks_2500-nodes_0.700000-pert';
    'BlockAcyclicPerturbation/8-blocks_2500-nodes_0.800000-pert';
    'BlockAcyclicPerturbation/8-blocks_2500-nodes_0.900000-pert';
    %'BlockAcyclicPerturbation/10-blocks_2500-nodes_0.000000-pert';
    %'BlockAcyclicPerturbation/10-blocks_2500-nodes_0.100000-pert';
    %'BlockAcyclicPerturbation/10-blocks_2500-nodes_0.200000-pert';
    %'BlockAcyclicPerturbation/10-blocks_2500-nodes_0.300000-pert';
    %'BlockAcyclicPerturbation/10-blocks_2500-nodes_0.400000-pert';
    %'BlockAcyclicPerturbation/10-blocks_2500-nodes_0.500000-pert';
    %'BlockAcyclicPerturbation/10-blocks_2500-nodes_0.600000-pert';
    %'BlockAcyclicPerturbation/10-blocks_2500-nodes_0.700000-pert';
    %'BlockAcyclicPerturbation/10-blocks_2500-nodes_0.800000-pert';
    %'BlockAcyclicPerturbation/10-blocks_2500-nodes_0.900000-pert';
};

num_trials       = 5;             % How many times to run BAS on each graph
beta             = 1.1;           % Parameter used for diffusion transition matrix
plotFlag         = false;         % Plot eigenvalues inside BAS
verbose          = false;         % BAS verbosity
estimate_blocks  = false;          % Whether to estimate k or fix it manually


%% Initialize Results
nc = length(graph_cases);
results = struct();

fprintf('========================================================\n');
fprintf('      Block-Acyclic Spectral Clustering (BAS) Benchmark \n');
fprintf('========================================================\n\n');


%% Loop over cases
for c = 1:nc
    fprintf('--------------------------------------------------------\n');
    fprintf('Processing graph file: %s\n', graph_cases{c});

    %% Step 1: Load graph from .mat file
    file_path     = sprintf("%s/Data/Synthetic/%s", pwd, graph_cases{c});
    file_path_adj = sprintf("%s_adj.mat", file_path);
    file_path_lab = sprintf("%s_lab.mat", file_path);

    W      = load(file_path_adj).W;
    labels = load(file_path_lab).labels;

    % Extract information
    info = split(graph_cases{c},"/");
    info = split(info(2),"_");

    k       = 0; % Number of blocks
    n       = 0; % Number of nodes
    epsilon = 0; % Perturbation magnitude

    for info_idx = 1:size(info,1)
        value = split(info(info_idx), "-");
        
        if strcmp(value(2), "blocks")
            k = str2double(value(1));
        elseif strcmp(value(2), "nodes")
            n = str2double(value(1));
        elseif strcmp(value(2), "pert")
            epsilon = str2double(value(1));
        end
    end

    n_edges = nnz(W);
    
    fprintf('Blocks                     : %d\n', k);
    fprintf('Nodes                      : %d\n', n);
    fprintf('Edges                      : %d\n', n_edges);
    fprintf('Magnitude of perturbation  : %f\n', epsilon);

    %% Initialize metrics storage for this graph
    RCut_all       = zeros(num_trials, 1);
    NCut_all       = zeros(num_trials, 1);
    NMI_all        = zeros(num_trials, 1);
    FScore_all     = zeros(num_trials, 1);
    Modularity_all = zeros(num_trials, 1);
    time_all       = zeros(num_trials, 1);

    RCut_half_all       = zeros(num_trials, 1);
    NCut_half_all       = zeros(num_trials, 1);
    NMI_half_all        = zeros(num_trials, 1);
    FScore_half_all     = zeros(num_trials, 1);
    Modularity_half_all = zeros(num_trials, 1);
    time_half_all       = zeros(num_trials, 1);

    RCut_diffusion_all       = zeros(num_trials, 1);
    NCut_diffusion_all       = zeros(num_trials, 1);
    NMI_diffusion_all        = zeros(num_trials, 1);
    FScore_diffusion_all     = zeros(num_trials, 1);
    Modularity_diffusion_all = zeros(num_trials, 1);
    time_diffusion_all       = zeros(num_trials, 1);

     %% Step 2: Run BAS clustering multiple times
    for t = 1:num_trials
        fprintf('Trial %d/%d...\n', t, num_trials);
        tStart = tic;

        %% Run BAS clustering
        [clusters, centroids] = BAS(W, k, "transition", 1, plotFlag, verbose, graph_cases{c});
        
        %% Compute partitioning metrics
        normalized = 1;
        [RCut, NCut, NMI, FScore, Modularity] = ComputeMetrics(labels, clusters, W);
        
        %% Store metrics
        RCut_all(t)       = RCut;
        NCut_all(t)       = NCut;
        NMI_all(t)        = NMI;
        FScore_all(t)     = FScore;
        Modularity_all(t) = Modularity;
        time_all(t)       = toc(tStart);

        tStart = tic;

        %% Run BAS clustering half
        [clusters, centroids] = BASnew(W, k, "transition", 1, plotFlag, verbose, graph_cases{c});
        
        %% Compute partitioning metrics
        normalized = 1;
        [RCut, NCut, NMI, FScore, Modularity] = ComputeMetrics(labels, clusters, W);
        
        %% Store metrics
        RCut_half_all(t)       = RCut;
        NCut_half_all(t)       = NCut;
        NMI_half_all(t)        = NMI;
        FScore_half_all(t)     = FScore;
        Modularity_half_all(t) = Modularity;
        time_half_all(t)       = toc(tStart);

        tStart = tic;

        %% Run BAS diffusion clustering
        [clusters, centroids] = BAS(W, k, "power", beta, plotFlag, verbose, graph_cases{c});
        
        %% Compute partitioning metrics
        normalized = 1;
        [RCut, NCut, NMI, FScore, Modularity] = ComputeMetrics(labels, clusters, W);
        
        %% Store metrics
        RCut_diffusion_all(t)       = RCut;
        NCut_diffusion_all(t)       = NCut;
        NMI_diffusion_all(t)        = NMI;
        FScore_diffusion_all(t)     = FScore;
        Modularity_diffusion_all(t) = Modularity;
        time_diffusion_all(t)       = toc(tStart);
    end

    %% Step 3: Save results for this graph
    results(c).graph_name           = graph_cases{c};
    results(c).perturbation         = epsilon;
    results(c).n_nodes              = n;
    results(c).n_edges              = n_edges;
    results(c).RCut                 = RCut_all;
    results(c).NCut                 = NCut_all;
    results(c).NMI                  = NMI_all;
    results(c).FScore               = FScore_all;
    results(c).Modularity           = Modularity_all;
    results(c).time                 = time_all;
    results(c).RCut_half            = RCut_half_all;
    results(c).NCut_half            = NCut_half_all;
    results(c).NMI_half             = NMI_half_all;
    results(c).FScore_half          = FScore_half_all;
    results(c).Modularity_half      = Modularity_half_all;
    results(c).time_half            = time_half_all;
    results(c).RCut_diffusion       = RCut_diffusion_all;
    results(c).NCut_diffusion       = NCut_diffusion_all;
    results(c).NMI_diffusion        = NMI_diffusion_all;
    results(c).FScore_diffusion     = FScore_diffusion_all;
    results(c).Modularity_diffusion = Modularity_diffusion_all;
    results(c).time_diffusion       = time_diffusion_all;


    %% Step 4: Print per-graph summary
    fprintf('Results for graph: %s\n', graph_cases{c});
    fprintf('--------------------------------------------------------\n');
    fprintf('Nodes        : %d\n', n);
    fprintf('Blocks       : %d\n', k);
    fprintf('Edges        : %d\n', n_edges);
    fprintf('Perturbation : %d\n\n', epsilon);
    fprintf('------ BAS ------\n');
    fprintf('F-Score (avg): %.4f | Std: %.4f\n', mean(FScore_all), std(FScore_all));
    fprintf('NCut (avg)   : %.4f | Std: %.4f\n', mean(NCut_all), std(NCut_all));
    fprintf('Modularity   : %.4f | Std: %.4f\n', mean(Modularity_all), std(Modularity_all));
    fprintf('Time (avg)   : %.4fs | Std: %.4fs\n\n', mean(time_all), std(time_all));
    fprintf('------ BAS half ------\n')
    fprintf('F-Score (avg): %.4f | Std: %.4f\n', mean(FScore_half_all), std(FScore_half_all));
    fprintf('NCut (avg)   : %.4f | Std: %.4f\n', mean(NCut_half_all), std(NCut_half_all));
    fprintf('Modularity   : %.4f | Std: %.4f\n', mean(Modularity_half_all), std(Modularity_half_all));
    fprintf('Time (avg)   : %.4fs | Std: %.4fs\n\n', mean(time_half_all), std(time_half_all));
    fprintf('------ BAS diffusion ------\n')
    fprintf('F-Score (avg): %.4f | Std: %.4f\n', mean(FScore_diffusion_all), std(FScore_diffusion_all));
    fprintf('NCut (avg)   : %.4f | Std: %.4f\n', mean(NCut_diffusion_all), std(NCut_diffusion_all));
    fprintf('Modularity   : %.4f | Std: %.4f\n', mean(Modularity_diffusion_all), std(Modularity_diffusion_all));
    fprintf('Time (avg)   : %.4fs | Std: %.4fs\n\n', mean(time_diffusion_all), std(time_diffusion_all));
    fprintf('--------------------------------------------------------\n\n');
end


%% Plotting
C=[ 0         0    1.0000
         0    0.4980         0
    1.0000    0.6000         0
    0.6353    0.0784    0.1843
     0.1490    0.8588    0.5059
    .0000    0.000         0  
    0.4000    0.2000    1.0000
    1.0000    0.0000    0.0000    
    ]; 

style={'-*','-x','-+','-.o','-.s','-.d','--<','-->'};
color={C(1,:),C(2,:),C(3,:),C(4,:)};

% Extract values
epsilons   = [results.perturbation];

NCut             = mean([results.NCut]);
FScore           = mean([results.FScore]);
NCut_half        = mean([results.NCut_half]);
FScore_half      = mean([results.FScore_half]);
NCut_diffusion   = mean([results.NCut_diffusion]);
FScore_diffusion = mean([results.FScore_diffusion]);

%% NCut
figure

plot(epsilons, NCut, style{1},'color',C(1,:),'LineWidth',2,'MarkerSize',8);hold on;
plot(epsilons, NCut_half, style{2},'color',C(2,:),'LineWidth',2,'MarkerSize',8);hold on; 
plot(epsilons, NCut_diffusion, style{3},'color',C(3,:),'LineWidth',2,'MarkerSize',8);hold on; 

ax = gca;
ax.YAxis(1).Color = 'k';

legend({...
    'BAS', 'BAS half', 'BAS diffusion'...
     },'interpreter','latex','location','northwest');
 
 
xlabel('Magnitude of perturbation $\epsilon$','interpreter','latex');


ylabel('Normalized Cut','interpreter','latex');

xlim([min(epsilons),max(epsilons)]);

set(gca,'fontsize',30);
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on','XTick',epsilons);

tightfig;

set(gcf,'units','points','position',[10 10 300 200]*1.9);

%% F-Score

figure

plot(epsilons, FScore, style{1},'color',C(1,:),'LineWidth',2,'MarkerSize',8);hold on;
plot(epsilons, FScore_half, style{2},'color',C(2,:),'LineWidth',2,'MarkerSize',8);hold on;
plot(epsilons, FScore_diffusion, style{3},'color',C(3,:),'LineWidth',2,'MarkerSize',8);hold on; 

ax = gca;
ax.YAxis(1).Color = 'k';

legend({...
    'BAS', 'BAS half', 'BAS diffusion'...
     },'interpreter','latex','location','northwest');
 
 
xlabel('Magnitude of perturbation $\epsilon$','interpreter','latex');

ylabel('F-Score','interpreter','latex');


xlim([min(epsilons),max(epsilons)]);

set(gca,'fontsize',30);
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on','XTick',epsilons);

tightfig;

set(gcf,'units','points','position',[10 10 300 200]*1.9);