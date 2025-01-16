clear; clc; close all;

rng(42);

% Benchmark dimensions
num_blocks    = 8;
num_nodes     = 1000;
conn_prob     = 0.7;
block_weights = [0.17; 0.2; 0.05; 0.2; 0.14; 0.04; 0.07; 0.13];
%pert_prob  = [0.0; 0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.65; 0.7; 0.75; 0.8];

pert_prob  = [0.0; 0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9];
num_tests  = size(pert_prob,1);

% Compute uniform membership probability
%block_weights = ones(num_blocks, 1) .* (1/num_blocks);


% Initialize metrics vectors
NCut   = zeros(num_tests, 2);
RCut   = zeros(num_tests, 2);
NMI     = zeros(num_tests, 2);
FScore = zeros(num_tests, 2);

% Generate the unperturbed Block-Cycle graph
[W, labels] = GenBlockCycle(num_blocks, num_nodes, conn_prob, block_weights);

for i = 1:size(pert_prob,1)
    % Generate Block-Cyclic Graph
    fprintf("Generating Block Cycle Graph with %d blocks, %d nodes," + ...
        " %f perturbation...\n", num_blocks, num_nodes, pert_prob(i,1));

    % Add perturbation to the Graph
    Pert = AddPerturbation(W,pert_prob(i));

    fprintf("Running Spectral Clustering...\n");
    % Get clusters with BCS
    [cluster_index, ~] = BCS(Pert, num_blocks, false, false);
    
    % Compute and save metrics
    fprintf("   Computing and saving metrics...\n");
    [RCut(i,1), NCut(i,1), NMI(i,1), FScore(i,1)] = ComputeMetrics( ...
        labels,cluster_index,Pert);
    
    fprintf("Running SVD...\n");
    % Get clusters with RAsNMF
    [cluster_index, ~] = SVD_clustering(Pert, num_blocks);

    % Compute and save metrics
    fprintf("   Computing and saving metrics...\n");
    [RCut(i,2), NCut(i,2), NMI(i,2), FScore(i,2)] = ComputeMetrics( ...
        labels,cluster_index,Pert);
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

%% NCut
figure


plot(pert_prob, NCut(:,1), style{1},'color',C(1,:),'LineWidth',2,'MarkerSize',8);hold on;
plot(pert_prob, NCut(:,2), style{2},'color',C(2,:),'LineWidth',2,'MarkerSize',8);hold on; 
%loglog(Lambda,nnz_c, style{3},'color',C(2,:),'LineWidth',2,'MarkerSize',8);hold on; 

ax = gca;
ax.YAxis(1).Color = 'k';



legend({...
    'BCS', 'SVD'...
     },'interpreter','latex','location','northwest');
 
 
xlabel('Perturbation Probability $\gamma$','interpreter','latex');


ylabel('Normalized Cut','interpreter','latex');


xlim([min(pert_prob),max(pert_prob)]);

set(gca,'fontsize',14);
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on','XTick',pert_prob);

tightfig;

set(gcf,'units','points','position',[10 10 300 200]*1.9);
saveas(gcf,'NCut_perturbation','pdf');


%% F-Score

figure

plot(pert_prob, FScore(:,1), style{1},'color',C(1,:),'LineWidth',2,'MarkerSize',8);hold on;
plot(pert_prob, FScore(:,2), style{2},'color',C(2,:),'LineWidth',2,'MarkerSize',8);hold on; 
%loglog(Lambda,nnz_c, style{3},'color',C(2,:),'LineWidth',2,'MarkerSize',8);hold on; 

ax = gca;
ax.YAxis(1).Color = 'k';

legend({...
    'BCS', 'SVD'...
     },'interpreter','latex','location','northwest');
 
 
xlabel('Perturbation Probability $\gamma$','interpreter','latex');

ylabel('F-Score','interpreter','latex');


xlim([min(pert_prob),max(pert_prob)]);

set(gca,'fontsize',14);
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on','XTick',pert_prob);

tightfig;

set(gcf,'units','points','position',[10 10 300 200]*1.9);
saveas(gcf,'FScore_perturbation','pdf');

