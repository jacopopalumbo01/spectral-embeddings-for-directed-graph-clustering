clear; clc; close all;
rng(47);

addpaths_SC(); % Make sure your BAS functions are in path

fprintf('========================================================\n');
fprintf('      Diffusion Transition Matrix Benchmark \n');
fprintf('========================================================\n\n');

%% Parameters
num_trials = 5;

best_value_evaluation = {
    'BlockAcyclicPerturbation/3-blocks_2500-nodes_0.400000-pert';
    'BlockAcyclicPerturbation/5-blocks_2500-nodes_0.400000-pert';
    'BlockAcyclicPerturbation/8-blocks_2500-nodes_0.400000-pert';
    'BlockAcyclicPerturbation/10-blocks_2500-nodes_0.400000-pert';
};

% Betas to evaluate
betas = (0:0.1:2)';
nc    = length(best_value_evaluation);
nb    = size(betas,1);

%% Cyclic case
fprintf("Evaluating block-acyclic graph\n");

% Initialize results
modularities = zeros(nc,nb);
fscores      = zeros(nc,nb);

for c = 1:nc
    fprintf('--------------------------------------------------------\n');
    fprintf('Processing graph file: %s\n', best_value_evaluation{c});

    %% Step 1: Load graph from .mat file
    file_path     = sprintf("%s/Data/Synthetic/%s", pwd, best_value_evaluation{c});
    file_path_adj = sprintf("%s_adj.mat", file_path);
    file_path_lab = sprintf("%s_lab.mat", file_path);

    W      = load(file_path_adj).W;
    labels = load(file_path_lab).labels;

    % Extract information
    info = split(best_value_evaluation{c},"/");
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

    % Initialize outputs
    modularity = zeros(nb, num_trials);
    fscore     = zeros(nb, num_trials);
    
    %% Step 2: Run BCS clustering multiple times with different betas
    for i = 1:nb
        fprintf("Evaluating beta=%f\n", betas(i));
        for t = 1:num_trials
            fprintf('Trial %d/%d...\n', t, num_trials);
            
            % Perform BCS using diffusion transition matrix
            [inferred_labels,~] = BAS(W, k, "power", betas(i));
            
            % Compute metrics
            [RCut, NCut, NMI, FScore, modularity] = ComputeMetrics(labels,inferred_labels,W);
            
            % Save temporary metrics
            modularity(i,t) = modularity;
            fscore(i,t)     = FScore;
        end
    end

    % Save metrics
    modularities(c,:) = mean(modularity,2);
    fscores(c,:)      = mean(fscore,2);
end

% Calculate averages
modularities_avg = mean(modularities);
fscores_avg      = mean(fscores);

% Find the best beta
[maxFscore, maxFscoreInd] = max(fscores_avg);
fprintf('----------------------------------------------------------\n');
fprintf("The optimal value is beta=%f, which leads to FScore=%f\n", betas(maxFscoreInd), maxFscore);
fprintf('----------------------------------------------------------\n');

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

%% Modularity
figure

plot(betas, modularities_avg, style{1},'color',C(1,:),'LineWidth',2,'MarkerSize',8);hold on;

ax = gca;
ax.YAxis(1).Color = 'k';

legend({...
    'Modularity'...
     },'interpreter','latex','location','northwest');
 
 
xlabel('$\beta$ factor - BAS','interpreter','latex');


ylabel('Modularity','interpreter','latex');

xlim([min(betas),max(betas)]);

set(gca,'fontsize',30);
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on','XTick',betas);

tightfig;

set(gcf,'units','points','position',[10 10 300 200]*1.9);

%% F-Score

figure

plot(betas, fscores_avg, style{1},'color',C(1,:),'LineWidth',2,'MarkerSize',8);hold on;

ax = gca;
ax.YAxis(1).Color = 'k';

legend({...
    'F-Score'...
     },'interpreter','latex','location','northwest');
 
 
xlabel('$\beta$ factor - BAS','interpreter','latex');

ylabel('F-Score','interpreter','latex');


xlim([min(betas),max(betas)]);

set(gca,'fontsize',30);
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on','XTick',betas);

tightfig;

set(gcf,'units','points','position',[10 10 300 200]*1.9);


%{
fprintf("Evaluating block-acyclic graph\n");
% Generate block connection probability
P = ConnectionProbabilityMatrix("acyclic", k, conn_prob);

% Generate graph
[W, nodes] = GenerateGraph(n,k,rho,P,0,epsilon);

num_experiments = size(betas,1);

modularities = zeros(num_experiments,1);

for i = 1:num_experiments
    [inferred_labels,~] = BAS(W, k, "power", betas(i));
    
    % Compute metrics
    [RCut, NCut, NMI, FScore, modularity] = ComputeMetrics(nodes, ...
        inferred_labels,W);
    modularities(i) = modularity;
    fscores(i)      = FScore;

    fprintf("   beta=%f -> modularity=%f FScore=%f\n", betas(i), ...
        modularities(i), fscores(i));
end

%% Modularity
figure

plot(betas, modularities, style{1},'color',C(1,:),'LineWidth',2,'MarkerSize',8);hold on;

ax = gca;
ax.YAxis(1).Color = 'k';

legend({...
    'Modularity'...
     },'interpreter','latex','location','northwest');
 
 
xlabel('$\beta$ factor - BAS','interpreter','latex');


ylabel('Modularity','interpreter','latex');

xlim([min(betas),max(betas)]);

set(gca,'fontsize',30);
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on','XTick',betas);

tightfig;

set(gcf,'units','points','position',[10 10 300 200]*1.9);

%% F-Score

figure

plot(betas, fscores, style{1},'color',C(1,:),'LineWidth',2,'MarkerSize',8);hold on;

ax = gca;
ax.YAxis(1).Color = 'k';

legend({...
    'F-Score'...
     },'interpreter','latex','location','northwest');
 
 
xlabel('$\beta$ factor - BAS','interpreter','latex');

ylabel('F-Score','interpreter','latex');


xlim([min(betas),max(betas)]);

set(gca,'fontsize',30);
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on','XTick',betas);

tightfig;

set(gcf,'units','points','position',[10 10 300 200]*1.9);
%}