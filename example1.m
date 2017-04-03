function [err] = example1( mode, dataset, n_samples, k_values , m, epsilon, algo_type, improvements)
% Applies insertion of a new (random) node from dataset to the rest of the
% dataset and compares methods to update the eigendecomposition of graph
% Laplacian.  If mode == 'plot' then also plots the dependency with k of
% svd(L1 - L0).
%
% input:
% mode - plots the graph iff mode == 'plot'
% dataset - name of the dataset wanted (MNIST, poker or yeast) or a matlab data matrix for user dataset
% n_samples - number of samples (rows of dataset) wanted
% k_values - vecor of k values for the graph Laplcian
% m - number of known eigenpairs
% epsilon - parameter for the Gaussian kernel
% algo_type: cell array. which alorhim to use for the update, see "error_and_time.m"
% improvements: see 'error_and_time.m'
% 
% output:
% err = size(algo_type) * size(k_values) matrix with columns beign algorithms used and rows the corresponding k-value of the graph Laplacian

    % parameters and variables
    eig_vals = zeros(4, size(k_values,2)); % matrix to store the 6 cacluted eignvalues of dL
    dims = zeros(size(k_values,2), 1); % array to store the size of dL
    eigvals_error = zeros(size(k_values,2), size(algo_type,2));
    eigvecs_error = zeros(size(k_values,2), size(algo_type,2));
    eigvals_corr_error = zeros(size(k_values,2), size(algo_type,2));
    eigvecs_corr_error = zeros(size(k_values,2), size(algo_type,2));

    err = zeros(size(algo_type,2), 4);
    
    % load dataset
    if strcmpi(dataset,'poker') == 1
        fprintf('dataset: poker\n');        
        X = csvread('poker-hand-training-true.data'); 
        n_samples = min(size(X,1), n_samples);
        X = X(randperm(size(X,1)), 1:(end - 1)); X = X(1:n_samples, :);
    elseif strcmpi(dataset,'MNIST') == 1
        fprintf('dataset: MNIST\n');
        X = loadMNISTImages('train-images.idx3-ubyte');X = X';X = X(randperm(60000),:);X=X(1:n_samples,:);
    elseif strcmpi(dataset,'yeast') == 1
        X = csvread('yeast.csv'); 
        X = X(randperm(size(X,1)),:);
        X = X(1:n_samples, 1:(end - 1));
    else
        fprintf('dataset: user\n');
        n_samples = min(size(dataset,1), n_samples);
        X = dataset(1:n_samples, :);
    end
    
    % check number of examples vs. k
    if max(k_values) > n_samples
        fprintf('max k value is bigger than number of samples provided. Aborting.\n'); 
        return;
    end
    
    % find kNN
    fprintf('performing kNN search... '); 
    IDX = knnsearch(X,X, 'K', max(k_values));
    fprintf('done.\n');    
    
    % build graph Laplacian
    
    for i=1:size(k_values,2)

        fprintf('k = %d\n', k_values(i));
        [L0,L1,min_k(i)] = build_graph_laplacian(X,  IDX(:,1:k_values(i)), epsilon);
        
        dL = L1 - L0;
       
        % removing empty rows\columns of dL to discover its true size
        dL( ~any(dL,2), : ) = [];
        dL( :, ~any(dL,1) ) = [];
        
        % saving info
        eig_vals(:,i) = svds(dL, 4); % eigenvalues
        dims(i) = size(dL,1); % size
        
        if strcmpi(mode, 'plot') == 0    
                [~,  err] = error_and_time({L0, L1}, m, 1,algo_type, improvements);
        end
    end
    
    

    if strcmpi(mode, 'plot') == 1    
            %%% plotting
            first = eig_vals(1,:);
            others = eig_vals(2:end,:);

            % log(1 - largest_eigenvalue) vs. log(k)
            figure;
            scatter(log2(min_k), log2(1 - first));
            xlabel('log(k)');
            ylabel('log(1 - largest singular value)');

            % log(eigenvalue) vs. log(k)
            figure;
            hold on;
            cmap = hsv(3);
            chars = ['o', '+', '*'];
            for j=1:size(others, 1)
                scatter(log2(min_k), log2(others(j,:)), 50, cmap(j,:), chars(j));
            end
            legend('2nd','3rd', '4th', '5th', '6th');
            xlabel('log(k)');
            ylabel('log(ith singular value)');
            hold off;

            % log(size(dL)) vs. log(k)
            figure;
            scatter(log2(min_k), log2(dims));
            xlabel('log(k)');
            ylabel('log(size(dL))');
            
    else
             return;
    
    end
    
end

