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
    scatter(real(cycle_eigval), imag(cycle_eigval), "ro");
    hold on;
    scatter(real(eigval), imag(eigval), "gx");
    legend('Cycle Eigenvalues', 'Other Eigenvalues');
    xlabel('Real part');
    ylabel('Imaginary part');
    title(eigval_title);

    % Plot cycle eigenvectors
    figure;
    for i = 1:size(cycle_eigvec, 2)
        plot(real(cycle_eigvec(:, i)), imag(cycle_eigvec(:, i)), 'r.', "MarkerSize", 18);
        hold on;
    end
    legend("Cycle Eigenvectors")
    xlabel('Real part');
    ylabel('Imaginary part');
    title(eigvec_title);
end
