clear;clc;
rng(47);

fprintf("--------------------------------------------\n");
fprintf("BENCHMARK FOR ESTIMATION OF NUMBER OF BLOCKS\n");
fprintf("--------------------------------------------\n");

%% Parameters
n         = 200;                            % Number of blocks
k         = [2; 3; 4; 5; 6; 7; 8; 9; 10;];  % Number of clusters
conn_prob = 0.7;                            % Connection probability between blocks

%%
num_tests = size(k,1);

silhouette_estimates = zeros(num_tests,1);
modularity_estimates = zeros(num_tests,1);
eigengaps_estimates  = zeros(num_tests,1);

% Loop over all cases
for i = 1:num_tests
    % Not uniform block membership distribution
    rho_not_uniform = rand(k(i),1);
    rho_not_uniform = rho_not_uniform / sum(rho_not_uniform);

    % Generate block connection probability
    P = ConnectionProbabilityMatrix("cyclic", k(i), conn_prob);
    
    % Generate graph
    [W, nodes] = StochasticBlockmodel(n,k(i),rho_not_uniform,P);

    % Compute estimations
    silhouette_estimates(i) = EstimateNumBlocksCyclic(W, n / 10);
    modularity_estimates(i) = EstimateNumBlocksCyclicWithModularity(W, n/10);
    eigengaps_estimates(i)  = EstimateNumBlocksCyclicWithEigengap(W, n/10);
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


figure;

plot(2:10, k,style{1},'color',C(1,:),'LineWidth',2,'MarkerSize',8);hold on;
plot(2:10, silhouette_estimates, style{2},'color',C(2,:),'LineWidth',2,'MarkerSize',8);hold on;
plot(2:10, modularity_estimates, style{3},'color',C(3,:),'LineWidth',2,'MarkerSize',8);hold on;
plot(2:10, eigengaps_estimates, style{4},'color',C(4,:),'LineWidth',2,'MarkerSize',8);hold on;

ax = gca;
ax.YAxis(1).Color = 'k';

legend({...
    'Ground Truth', 'Silhouette', 'Modularity', 'Eigengaps'...
     },'interpreter','latex','location','northwest');
 
 
xlabel('Real Number','interpreter','latex');


ylabel('Estimated Number','interpreter','latex');

xlim([2,10]);

set(gca,'fontsize',20);
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on','XTick',2:10);

tightfig;

set(gcf,'units','points','position',[10 10 300 200]*1.9);
