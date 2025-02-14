clear;clc;
rng(47);

fprintf("-------------------------------------------\n");
fprintf("BENCHMARK FOR BAS ON STOCHASTIC BLOCKMODELS\n\n");
fprintf("-------------------------------------------\n");

%% Parameters
n   = 1000;         % Number of nodes
k   = 8;            % Number of blocks  
P1   = zeros(k,k);  % Block connection probability

% Block membership distribution
rho = [0.18; 0.2; 0.05; 0.2; 0.14; 0.04; 0.07; 0.13];
%rho = ones(k,1);
%rho = rho / sum(rho);
% Set block connection probability
conn_prob = 0.5;
for i = 1:k
    for j = 1:k
        if i < j
            P1(i,j) = conn_prob;
        end
    end
end

%% Generation of unperturbed graph
fprintf("Generating unperturbed graph\n")
[A1, tau1] = StochasticBlockmodel(n,k,rho,P1);


%% Testing phase

% Perturbation probabilities
epsilons  = [0.0; 0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9];
num_tests  = size(epsilons,1);

% Initialize metrics vectors
NCut   = zeros(num_tests, 2);
RCut   = zeros(num_tests, 2);
NMI     = zeros(num_tests, 2);
FScore = zeros(num_tests, 2);

% Create perturbing block connection probability matrix
Q  = rand(k,k);

for i = 1:num_tests
    fprintf("Testing epsilon=%f\n", epsilons(i,1));

    % Computing new perturbing block connection probability matrix given
    % the magnitude of perturbation epsilon
    P2 = epsilons(i,1) * Q;

    % Generate perturbing graph
    fprintf("   Generating perturbing graph\n");
    [A2, ~] = StochasticBlockmodel(n,k,rho,P2);

    % Combine unperturbed and perturbing graph
    fprintf("   Combining unperturbed and perturbing graph\n");
    A = CombineBlockmodels(A1,A2);

    fprintf("   Running Spectral Clustering\n");
    % Get clusters with BCS
    [cluster_index, ~] = BAS(A, k, false, false);
    
    % Compute and save metrics
    fprintf("       Computing and saving metrics\n");
    [RCut(i,1), NCut(i,1), NMI(i,1), FScore(i,1)] = ComputeMetrics( ...
        tau1, cluster_index, A);
    
    fprintf("   Running SVD\n");
    % Get clusters with RAsNMF
    [cluster_index, ~] = SVD_clustering(A, k);

    % Compute and save metrics
    fprintf("       Computing and saving metrics\n");
    [RCut(i,2), NCut(i,2), NMI(i,2), FScore(i,2)] = ComputeMetrics( ...
        tau1,cluster_index,A);

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

plot(epsilons, NCut(:,1), style{1},'color',C(1,:),'LineWidth',2,'MarkerSize',8);hold on;
plot(epsilons, NCut(:,2), style{2},'color',C(2,:),'LineWidth',2,'MarkerSize',8);hold on; 
%loglog(Lambda,nnz_c, style{3},'color',C(2,:),'LineWidth',2,'MarkerSize',8);hold on; 

ax = gca;
ax.YAxis(1).Color = 'k';

legend({...
    'BAS', 'SVD'...
     },'interpreter','latex','location','northwest');
 
 
xlabel('Magnitude of perturbation $\epsilon$','interpreter','latex');


ylabel('Normalized Cut','interpreter','latex');

xlim([min(epsilons),max(epsilons)]);

set(gca,'fontsize',50);
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on','XTick',epsilons);

tightfig;

set(gcf,'units','points','position',[10 10 300 200]*1.9);
saveas(gcf,'NCut_perturbation','pdf');

%% F-Score

figure

plot(epsilons, FScore(:,1), style{1},'color',C(1,:),'LineWidth',2,'MarkerSize',8);hold on;
plot(epsilons, FScore(:,2), style{2},'color',C(2,:),'LineWidth',2,'MarkerSize',8);hold on; 
%loglog(Lambda,nnz_c, style{3},'color',C(2,:),'LineWidth',2,'MarkerSize',8);hold on; 

ax = gca;
ax.YAxis(1).Color = 'k';

legend({...
    'BAS', 'SVD'...
     },'interpreter','latex','location','northwest');
 
 
xlabel('Magnitude of perturbation $\epsilon$','interpreter','latex');

ylabel('F-Score','interpreter','latex');


xlim([min(epsilons),max(epsilons)]);

set(gca,'fontsize',50);
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on','XTick',epsilons);

tightfig;

set(gcf,'units','points','position',[10 10 300 200]*1.9);
saveas(gcf,'FScore_perturbation','pdf');



