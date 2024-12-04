function PlotCyclicEig(eigval, cycle_eigval, cycle_eigvec)
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
    %
    figure;
    scatter(real(cycle_eigval), imag(cycle_eigval), "ro");
    hold on;
    scatter(real(eigval), imag(eigval), "gx");
    legend('Cycle Eigenvalues', 'Other Eigenvalues');
    xlabel('Real part');
    ylabel('Imaginary part');
    title("Eigenvalues of Transition Matrix");

    % Plot cycle eigenvectors
    figure;
    for i = 1:size(cycle_eigvec, 2)
        plot(real(cycle_eigvec(:, i)), imag(cycle_eigvec(:, i)), 'r.', "MarkerSize", 18);
        hold on;
    end
    legend("Cycle Eigenvectors")
    xlabel('Real part');
    ylabel('Imaginary part');
    title("Cycle eigenvectors");
end
