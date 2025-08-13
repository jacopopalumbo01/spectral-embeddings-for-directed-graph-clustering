% Path to results
results_path = "experiments/results/DSBM_increasing_clusters.mat";

% Select metrics you want to plot
metrics = {
    "Fscore",
};

% Select methods you want to plot
methods = {
     "SVD_unscaled",
     "SVD_unscaled_tSNE",
    % "SVD_scaled",
    % "SVD_scaled_tSNE",
    % "SKEW",
    % "SKEW_tSNE",
    % "Herm",
    % "Herm_tSNE",
    % "BAS",
    % "BAS_tSNE",
    % "BCS",
    % "BCS_tSNE"
};

% Select steps you want to plot
steps = {
    "DSBM_2blocks",
    "DSBM_3blocks",
    "DSBM_4blocks",
    "DSBM_5blocks",
    "DSBM_6blocks",
    "DSBM_7blocks",
    "DSBM_8blocks",
    "DSBM_9blocks",
    "DSBM_10blocks",
};

% Number of steps displayed
num_steps = [2;3;4;5;6;7;8;9;10];

%% Load matrix
res = load(results_path).results;

%% Extract needed values
metrics_values = {};

% for each metric
for i = 1:length(metrics)
    % for each graph
    for j = 1:length(steps)
        % for each method
        for k = 1:length(methods)
            % get strings
            step = string(steps(j));
            method = string(methods(k));
            metric = string(metrics(i));

            % compute average
            metrics_values.(metric).(method).avg(j) = mean(res.(step).(method).(metric));
            
            % compute minimum
            metrics_values.(metric).(method).min(j) = min(res.(step).(method).(metric));   

            % compute maximum
            metrics_values.(metric).(method).max(j) = max(res.(step).(method).(metric));
        end
    end
end


%% Plotting
C=[ 0         0.5     1.0000
    0         0       1.0000
    1.0000    0.6000  0
    1.0000    0.2000  0
  ]; 

style={'--*','-*','--x','-x'};

% for each metric
for i = 1:length(metrics)
    metric = string(metrics(i));

    % Start plotting
    figure;
    save_path = "DSBM";

    % for each method
    for j = 1:length(methods)
        method = string(methods(j));
        mvals  = metrics_values.(metric).(method);
        save_path = sprintf("%s-%s", save_path, method);

        % Plot confidence interval
        ylow = mvals.min;
        yhigh = mvals.max;

        % Shaded band (proper polygon)
        xconf = [num_steps(:); flipud(num_steps(:))];
        yconf = [ylow(:);      flipud(yhigh(:))];
        p = fill(xconf, yconf, C(j,:), 'HandleVisibility','off');
        p.FaceAlpha = 0.1;
        p.EdgeColor = 'none';  
        hold on;

        % Plot average line
        plot(num_steps, metrics_values.(metric).(method).avg, style{j},'color',C(j,:),'LineWidth',2,'MarkerSize',8,'DisplayName',strrep(method,"_"," "));
        hold on;
    end
    % Set legend
    legend("Interpreter","latex",'location','northeast');
    legend boxoff;

    % Set labels
    xlabel('Number of clusters','interpreter','latex');
    ylabel(metric,'interpreter','latex');

    xlim([min(num_steps),max(num_steps)]);

    set(gca,'fontsize',30);
    set(gca,'YMinorTick','on')
    set(gca,'XMinorTick','on','XTick',num_steps);
    
    tightfig;
end