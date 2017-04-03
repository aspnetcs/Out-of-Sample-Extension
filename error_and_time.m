function [timing, err,  algo_vals, new_vals] = error_and_time(matrix_type,  m, n_rounds, algo_type, improvements)
 
% returns the error and runtime of n_rounds independent experiments for
% rank-one update.
% matrix_type : if number, randomizes a matrix and a vector of size
% matrix_type to do the updates on. Otherwise, should be of the form
% {L0,L1} when L0 is a graph Laplacian before insertion of a new node and
% L1 is the same graph Laplacian after the insertion. 
% m = number of known eigenparis
% algo_type = a cell array consisting of string: 
%       'no_update' - use the old eigendeomposition to aproximate the new
%       'nystrom' - use nystrom method
%       'full' - use suggested algorithm as if all eigendecomosition is known
%        '10'/'11'/'12'/'20'/'21'/'22' - see 'update_eigenspectrum.m'
% improvements = a cell array 'pertrubation_correction'} if you wish to
% apply peturbation correction to any of the above methods.
          
    eig_vals_error = zeros(n_rounds, size(algo_type,2));
    eig_vecs_error = zeros(n_rounds, size(algo_type,2));
    eig_vals_corr_error = zeros(n_rounds, size(algo_type,2));
    eig_vecs_corr_error = zeros(n_rounds, size(algo_type,2));
    
    timing = zeros(n_rounds, size(algo_type,1) + 1); % 1st columns matlab, 2+ column algo

    for round = 1:n_rounds
       
        is_GL = ~(size(matrix_type,1) == 1 && size(matrix_type,2) == 1); % 1 <=> graph Laplcian
        
        if is_GL == 0 % random matrix and vector
            fprintf('Using random matrix and vector.\n');
            n = matrix_type;
            A = randn(n); A = A'*A; A = A/norm(A); 
            [Q,D] = eig(A); D = diag(D); [~,I] = sort(D); D(I(1:(n-m))) = -10  + rand(n-m,1)*1e-12; A = Q*real(diag(D))*Q';
            v = 1 + randn(n, 1); % I = randperm(n);  I = I(1:n - 100); v(I) = 0;
            v = v/norm(v);
            lambda = 1;
            B = A  + lambda * (v * v');
            [org_vecs, org_vals] = eigs(A + 10*eye(n), m);
            [~,P] = sort(diag(org_vals), 'descend'); org_vecs = org_vecs(:,P); org_vals = org_vals(P,P);      


        else %dataset matrix and Graph Laplacian
            A = matrix_type{1};
            n = size(A,1);
            B = matrix_type{2};
            [v, lambda] = eigs(B - A, 1);
            [org_vecs,org_vals] = eigs(A(2:end,2:end) + 10*eye(n - 1), m - 1) ; [~,I] = sort(diag(org_vals), 'descend'); 
            org_vals = org_vals(I,I); org_vecs = org_vecs(:, I);
            org_vecs = [ 1, zeros(1,m - 1) ; zeros(n - 1, 1) , org_vecs];
            org_vals = [11, zeros(1,m - 1) ; zeros(m - 1, 1), org_vals];
            
        end
  

         % find the real decompostion for comparison
         tic; [new_vecs, new_vals] = eigs(B + 10*eye(n), m); timing(round, 1) = toc;
         [~,P] = sort(diag(new_vals), 'descend'); new_vecs = new_vecs(:,P); new_vals = new_vals(P,P);
        C = B - (A + lambda * (v * v'));
        
        if is_GL, new_vals = new_vals(2:end-1,2:end-1); new_vecs = new_vecs(:, 2:end-1); end
        
        for i = 1:size(algo_type,2)
            type = char(algo_type(i));
            fprintf('[+] Using method: %s\n', type);
            if strcmpi(type, 'no_update')
                algo_vals = diag(org_vals);
                algo_vecs = org_vecs;
                if is_GL, algo_vals = algo_vals(3:end); algo_vecs = algo_vecs(:, 3:end); end
                 
            elseif strcmpi(type, 'nystrom')
                tic;
                algo_vals = (n/(n-1)) * (diag(org_vals)-10) + 10;
                for j = 1:m
                    algo_vecs(:, j) =sqrt((n-1)/n) * (1/(org_vals(j, j) - 10)) * B(:, 2:end) * org_vecs(2:end,j);
                end
                timing(round, i + 1) = toc;
                if is_GL, algo_vals = algo_vals(3:end); algo_vecs = algo_vecs(:, 3:end); end
                
            elseif strcmpi(type, 'full')
                tic;
                [algo_vecs,  algo_vals] = eigs(A + lambda * (v*v') + 10*eye(n), m); 
                [~,P] = sort(diag(algo_vals), 'descend'); algo_vecs = algo_vecs(:,P); algo_vals = algo_vals(P,P);
                algo_vals = diag(algo_vals);   
                timing(round, i + 1) = toc;
                if is_GL, algo_vals = algo_vals(2:end-1); algo_vecs = algo_vecs(:, 2:end-1); end
            else % our algo
                [algo_vecs, algo_vals, t] = update_eigenspectrum(A + 10*eye(n), lambda, v, org_vecs, org_vals , 0, type); 
                timing(round, i + 1) = t;
                if is_GL, algo_vals = algo_vals(3:end); algo_vecs = algo_vecs(:, 3:end); end
            end
            

            
          eig_vals_error(round, i) = mean(abs(algo_vals - diag(new_vals))./abs(diag(new_vals)));
          eig_vecs_error(round, i) = mean(min(sum((new_vecs - algo_vecs).^2,1), sum((new_vecs+ algo_vecs).^2,1)))^0.5;
            
           for k = 1:size(improvements,2)
               improvement = improvements(k);
              if strcmpi(improvement, 'pertrubation_correction') == 1
                            fprintf('Applying pertrubation correction... ');
                            [algo_vec_cor, algo_vals_cor] = pertrubation_correction(algo_vecs, algo_vals, C);
                            fprintf('Done. \n');
                            eig_vals_corr_error(round, i) = mean(abs(algo_vals_cor - diag(new_vals))./abs(diag(new_vals)));
                            eig_vecs_corr_error(round, i) = mean(min(sum((new_vecs - algo_vec_cor).^2,1), sum((new_vecs+ algo_vec_cor).^2,1)))^0.5;       
              end
           end
          
        end
        
        
        
    end

     eig_vals_error = mean(eig_vals_error, 1);
     eig_vecs_error = mean(eig_vecs_error, 1);
     eig_vals_corr_error = mean(eig_vals_corr_error, 1);
     eig_vecs_corr_error = mean(eig_vecs_corr_error, 1);
     
     err = [eig_vals_error' ,  eig_vals_corr_error', eig_vecs_error' ,eig_vecs_corr_error'];
     timing = mean(timing,1);
end

