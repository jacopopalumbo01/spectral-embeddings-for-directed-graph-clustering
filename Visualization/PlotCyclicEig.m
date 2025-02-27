function PlotCyclicEig(eigval, cycle_eigval, cycle_eigvec, graph_name)
% PlotCyclic - Plots block-cyclic eigenvalues and eigenvectors
%
%% Syntax:
%        PlotCyclicEig(eigval, cycle_eigval, cycle_eigvec)
%
%% Input Arguments:
%       *Required Input Arguments*
%       - eigval:           Other eigenvalues
%       - cycle_eigval:     Cycle eigenvalues
%       - cycle_eigvec:     Cycle eigenvectors
%       - graph_name:       Graph name
%

if nargin < 4
    eigval_title = "Eigenvalues of Transition Matrix";
    eigvec_title = "Cycle eigenvectors";
else
    eigval_title = strcat("Eigenvalues of Transition Matrix - ", graph_name);
    eigvec_title = strcat("Cycle eigenvectors - ", graph_name);
end

figure;
scatter(real(cycle_eigval), imag(cycle_eigval), 300, [0.4660 0.6740 0.1880], "o", "LineWidth",3);
hold on;
scatter(real(eigval), imag(eigval), 300, [0.8500 0.3250 0.0980], "x", "LineWidth",3);
legend('Cycle Eigenvalues', 'Other Eigenvalues', 'FontSize', 15);
xlabel('Real part', 'FontSize',15);
ylabel('Imaginary part', 'FontSize', 15);
title(eigval_title);

set(gca,'fontsize',15);

tightfig;


% Plot cycle eigenvectors
figure;
for i = 1:size(cycle_eigvec, 2)
    plot(real(cycle_eigvec(:, i)), imag(cycle_eigvec(:, i)), 'r.', "MarkerSize", 50, "LineWidth", 3);
    hold on;
end
legend("Cycle Eigenvectors", "FontSize", 15)
xlabel('Real part', "FontSize", 15);
ylabel('Imaginary part', "FontSize", 15);
title(eigvec_title);
set(gca,'fontsize',15);

tightfig;
end
