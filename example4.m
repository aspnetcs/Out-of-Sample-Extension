function [ err ] = example4( n_rounds, dataset, n_samples, k_values , m, epsilon, algo_type, improvements )
% return the mean error for several independent Laplcian eigendecompositon
% using several algorithms
% updates
% input:
% n_rounds = number of experiments to perform 
% dataset - name of the dataset wanted (MNIST, poker or yeast) or a matlab data matrix for user dataset
% n_samples - number of samples (rows of dataset) wanted
% k_values - vecor of k values for the graph Laplcian
% m - number of known eigenpairs
% epsilon - parameter for the Gaussian kernel
% algo_type: cell array. which alorhim to use for the update, see "error_and_time.m"
% improvements: see 'error_and_time.m'
% output:
% err = a size(algo_type)  * 4 matrix of the mean relative error. The
% columns correspond to the error of: eigenvalues, eigenvalues with
% correction, eigenvectors, eigenvectors with correction

        err = zeros(size(algo_type,2), 4);

        for round = 1:n_rounds
            fprintf('round %d/%d \n', round, n_rounds);
            err = err + example1('no_plot', dataset, n_samples, k_values , m, epsilon, algo_type, improvements);
        end
        
        err = err / n_rounds;
end

