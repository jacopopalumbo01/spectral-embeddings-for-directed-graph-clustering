clear;close all;

rng(47);

fprintf("------------------------------------------\n");
fprintf("BENCHMARK FOR ADJUSTABLE TRANSITION MATRIX\n");
fprintf("------------------------------------------\n");

%% Parameters
n         = 1000;         % Number of nodes
k         = 8;            % Number of blocks  
conn_prob = 0.8;          % Connection probability between blocks
epsilon   = 0.4;          % Magnitude of perturbation

% Block membership distribution
rho = [0.18; 0.2; 0.05; 0.2; 0.14; 0.04; 0.07; 0.13];

% Betas to evaluate
betas = (0:0.1:2)';

%% Cyclic case
fprintf("Evaluating block-cyclic graph\n");
% Generate block connection probability
P = ConnectionProbabilityMatrix("cyclic", k, conn_prob);

% Generate graph
[W, nodes] = GenerateGraph(n,k,rho,P,epsilon);

num_experiments = size(betas,1);

modularities = zeros(num_experiments,1);
fscores      = zeros(num_experiments,1);

for i = 1:num_experiments
    [inferred_labels,~] = BCS(W, k, "power", betas(i));
    
    % Compute metrics
    [RCut, NCut, NMI, FScore, modularity] = ComputeMetrics(nodes, ...
        inferred_labels,W);
    modularities(i) = modularity;
    fscores(i)      = FScore;

    fprintf("   beta=%f -> modularity=%f FScore=%f\n", betas(i), ...
        modularities(i), fscores(i));
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

%% Modularity
figure

plot(betas, modularities, style{1},'color',C(1,:),'LineWidth',2,'MarkerSize',8);hold on;

ax = gca;
ax.YAxis(1).Color = 'k';

legend({...
    'Modularity'...
     },'interpreter','latex','location','northwest');
 
 
xlabel('$\beta$ factor - BCS','interpreter','latex');


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
 
 
xlabel('$\beta$ factor - BCS','interpreter','latex');

ylabel('F-Score','interpreter','latex');


xlim([min(betas),max(betas)]);

set(gca,'fontsize',30);
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on','XTick',betas);

tightfig;

set(gcf,'units','points','position',[10 10 300 200]*1.9);



fprintf("Evaluating block-acyclic graph\n");
% Generate block connection probability
P = ConnectionProbabilityMatrix("acyclic", k, conn_prob);

% Generate graph
[W, nodes] = GenerateGraph(n,k,rho,P,epsilon);

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
