function [eig_vals] = example1( dataset, n_samples, k_values )

    % parameters and variables
    eig_vals = zeros(6, size(k_values,2)); % matrix to store the 6 cacluted eignvalues of dL
    dims = zeros(size(k_values,2), 1); % array to store the size of dL

    % load dataset
    if strcmpi(dataset,'poker') == 1
        fprintf('dataset: poker\n');        
        X = csvread('poker-hand-training-true.data'); 
        n_samples = min(size(X,1), n_samples);
        X = X(randperm(size(X,1)), 1:(end - 1)); X = X(1:n_samples, :);
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

        fprintf('calculating graph Laplacian matrices... '); 
        W1 = zeros(n_samples,n_samples);
        for m = 1:n_samples
            for n = 1:k_values(i)
                 curr = IDX(m,n);
                 W1(m,curr) = 1;
                 W1(curr,m) = 1;
            end
        end
        W1 = sparse(W1);
        W0 = W1; W0(1,:) = 0; W0(:,1) = 0; W0(1,1) = 1;
        Dinv0 = sparse(diag(1./sum(W0,2))); L0  = (Dinv0.^(0.5))*(W0)* (Dinv0.^(0.5));
        clear W0;clear Dinv0;
        D1 = sparse(diag(sum(W1,2))); Dinv1 = sparse(diag(1./sum(W1,2))); L1  = (Dinv1.^(0.5))*(W1)* (Dinv1.^(0.5));
        min_k(i) = min(diag(D1(2:end, 2:end))); 
        clear W1; clear D1; clear Dinv1;
        
        fprintf('done.\n');  

        dL = L1 - L0;
        
        % removing empty rows\columns of dL to discover its true size
        dL( ~any(dL,2), : ) = [];
        dL( :, ~any(dL,1) ) = [];

        % saving info
        eig_vals(:,i) = svds(dL); % eigenvalues
        dims(i) = size(dL,1); % size

    end
    
    %%% plotting
    first = eig_vals(1,:);
    others = eig_vals(2:end,:);
    
    % log(size(dL)) vs. log(k)
    figure;
    scatter(log10(min_k), log10(dims));
    xlabel('log(k)');
    ylabel('log(size(dL))');

    % log(1 - largest_eigenvalue) vs. log(k)
    figure;
    scatter(log10(min_k), log10(1 - first));
    xlabel('log(k)');
    ylabel('log(1 - largest singular value)');


    % log(eigenvalue) vs. log(k)
    figure;
    hold on;
    cmap = hsv(6);
    for j=1:size(others, 1)
        scatter(log10(min_k), log10(others(j,:)), 50, cmap(j,:), 'filled');
    end
    legend('2nd','3rd', '4th', '5th', '6th');
    xlabel('log(k)');
    ylabel('log(ith singular value)');
    hold off;
end

