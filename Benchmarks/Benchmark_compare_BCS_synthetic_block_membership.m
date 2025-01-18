clear; clc; close all;

rng(42);

% Benchmark dimensions
num_blocks    = 10;
num_nodes     = 1000;
conn_prob     = 0.7;
pert_prob  = 0.0;


% Compute Block Weights probabilities with different epsilon
mu = num_blocks;        % Average for our sampling. When epsilon = 0 
                        % we obtain a uniform distribution

% Initialize block weights
block_weights = zeros(num_blocks,10);

for i = 0:9 
    epsilon = i * 0.1;  % epsilon from 0 to 0.9
    
    % Create vector of positions from 0 to 1
    positions = linspace(0, 1, num_blocks)';
    
    % Create weights that are always positive by using positions from 0 to 1
    weights = ones(num_blocks, 1) + (epsilon * num_blocks * positions);
    
    % Store in block_weights
    block_weights(:,i + 1) = weights;
    
    % Normalize to sum = 1
    block_weights(:,i + 1) = block_weights(:,i + 1) / sum(block_weights(:,i + 1));
    
    % Compute MSE
    err = mean(((ones(num_blocks,1) ./ num_blocks) - block_weights(:,i + 1)).^2);
    
    fprintf("Epsilon=%f, MSE=%f\n", epsilon, err);
end


num_tests = 10;

% Initialize metrics vectors
NCut   = zeros(num_tests, 2);
RCut   = zeros(num_tests, 2);
NMI     = zeros(num_tests, 2);
FScore = zeros(num_tests, 2);


for i = 1:num_tests
    % Generate Block-Cyclic Graph
    fprintf("Generating Block Cycle Graph with %d blocks, %d nodes\n" ...
        , num_blocks, num_nodes);

    % Generate the unperturbed Block-Cycle graph
    [W, labels] = GenBlockCycle_dms(num_blocks, num_nodes, conn_prob, block_weights(:,i));


    fprintf("Running Spectral Clustering...\n");
    % Get clusters with BCS
    [cluster_index, ~] = BCS(W, num_blocks, false, false);
    
    % Compute and save metrics
    fprintf("   Computing and saving metrics...\n");
    [RCut(i,1), NCut(i,1), NMI(i,1), FScore(i,1)] = ComputeMetrics( ...
        labels,cluster_index,W);

    fprintf("Running SVD...\n");
    % Get clusters with RAsNMF
    [cluster_index, ~] = SVD_clustering(W, num_blocks);

    % Compute and save metrics
    fprintf("   Computing and saving metrics...\n");
    [RCut(i,2), NCut(i,2), NMI(i,2), FScore(i,2)] = ComputeMetrics( ...
        labels,cluster_index,W);
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

epsilon = transpose(linspace(0,9,10)) * 0.1;
plot(epsilon, NCut(:,1), style{1},'color',C(1,:),'LineWidth',2,'MarkerSize',8);hold on;
plot(epsilon, NCut(:,2), style{2},'color',C(2,:),'LineWidth',2,'MarkerSize',8);hold on; 
%loglog(Lambda,nnz_c, style{3},'color',C(2,:),'LineWidth',2,'MarkerSize',8);hold on; 

ax = gca;
ax.YAxis(1).Color = 'k';

legend({...
    'BCS', 'SVD'...
     },'interpreter','latex','location','northwest');
 
 
xlabel('Perturbation $\epsilon$','interpreter','latex');


ylabel('Normalized Cut','interpreter','latex');

xlim([0,1]);

set(gca,'fontsize',14);
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on','XTick',epsilon);

tightfig;

set(gcf,'units','points','position',[10 10 300 200]*1.9);
saveas(gcf,'NCut_perturbation','pdf');

pause
%% F-Score

figure

plot(epsilon, FScore(:,1), style{1},'color',C(1,:),'LineWidth',2,'MarkerSize',8);hold on;
plot(epsilon, FScore(:,2), style{2},'color',C(2,:),'LineWidth',2,'MarkerSize',8);hold on; 
%loglog(Lambda,nnz_c, style{3},'color',C(2,:),'LineWidth',2,'MarkerSize',8);hold on; 

ax = gca;
ax.YAxis(1).Color = 'k';

legend({...
    'BCS', 'SVD'...
     },'interpreter','latex','location','northwest');
 
 
xlabel('Perturbation $\epsilon$','interpreter','latex');

ylabel('F-Score','interpreter','latex');


xlim([0,1]);

set(gca,'fontsize',14);
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on','XTick',epsilon);

tightfig;

set(gcf,'units','points','position',[10 10 300 200]*1.9);
saveas(gcf,'FScore_perturbation','pdf');