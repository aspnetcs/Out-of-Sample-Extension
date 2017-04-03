function [algo_vecs, algo_vals, time] = update_eigenspectrum(A, lambda, v, org_vecs, org_vals, sparse, algo_type)
%	Update eigenspectrum of symmetric matrix A after rank-one
%	pertrubation, i.e., find the eigenspectrum of A_ = A + lambda * v'v;
%	Using the method described in the paper. If the original
%	eigendecompositin is not provided (org_vecs/org_vals) It will be
%	calulated.
% algo_type = a 2-char string of the form (0/1/2)(1/2) indicating the type of
% equation to solve. Left char indicates TSE - 1 or CTSE - 2. Right char
% indicates TEF - 0, CTEF1 - 1, CEEF2 - 2

    if strcmpi(algo_type, '11') == 1
        se_type = 'TSE';
        vf_type = 'TEF';
    elseif strcmpi(algo_type, '12') == 1
        se_type = 'TSE';
        vf_type = 'CTEF';
    elseif strcmpi(algo_type, '21') == 1
        se_type = 'CTSE';
        vf_type = 'TEF';
    elseif strcmpi(algo_type, '22') == 1
        se_type = 'CTSE';
        vf_type = 'CTEF';
    elseif strcmpi(algo_type, '10') == 1
        se_type = 'TSE';
        vf_type = 'none';
    elseif strcmpi(algo_type, '20') == 1
        se_type = 'CTSE';
        vf_type = 'none';
    else
        fprintf('unkown algorithm type. Aborting. \n');
        return;
    end
    
    % remove very small entries of v and normalize
    v(abs(v) < 1e-12) = 0;
    v = v/norm(v);
    
%     if lambda <= 0
%         fprintf('lambda must be positive. Aborting.');
%         return;
%     end
    
    if size(org_vecs, 2) ~= size(org_vals,1)
        fprintf('number of input eigevalues should be the same as number of input eigenvectors. Aborting.');
        return;
    end
    
    if sparse == 0 && issparse(A) == 1
        A = full(A);
        v = full(v);
    end
    
    if sparse == 1 && issparse(A) == 0
        A = sparse(A);
        v = sparse(v);
    end
    
    if size(org_vecs, 1) == 1 && size(org_vecs, 2) == 1 % org_eigvecs size tell how many pairs to calculate
        fprintf('old eigenvectors/old eigenvalues were not provided. Calculating the first %d eigenpairs... ', org_vecs);
        [org_vecs, org_vals] = eigs(A, org_vecs);
        fprintf('done.\n');
    end

    % descending order
    [~,P] = sort(diag(org_vals), 'descend'); org_vecs = org_vecs(:,P); org_vals = org_vals(P,P);
    
    % variables and parameters

    n = size(org_vecs, 1); % dimension
    m = size(org_vecs, 2); % number of known eigenpairs
    
    dd = diag(org_vals);
    algo_vals = zeros(m, 1);
    algo_vecs = zeros(n, m);
    prox = 1e-14;
    
   tic;
    z = org_vecs'*v;
    
   factor1 = ((1 - norm(z)^2));
   factor2 = ((1 - norm(z)^2)*(trace(A) - sum(dd)))/(n-m);
    % calculate eigenvalues
    if strcmpi(se_type,'CTSE') == 1 && m ~= n
        f = inline('1 + lambda * sum(sum((u.^2)./(d - x))) + lambda * (1 - sum(u.^2))./ (-x) + lambda * factor2./(-x^2)', 'x', 'u', 'd', 'factor2', 'lambda', 'factor1');
        df = inline('lambda * sum(sum((u.^2)./(d - x).^2)) + lambda * (1 - sum(u.^2))./ (x^2) + 2 * lambda *factor2./(x^3)', 'x', 'u', 'd', 'factor2', 'lambda', 'factor1');
    else
        f =  inline('1 + lambda * sum(sum((u.^2)./(d - x))) + lambda * factor1./ (-x) ', 'x', 'u', 'd', 'factor2', 'lambda', 'factor1');
        df = inline('lambda * sum(sum((u.^2)./((d - x).^2))) + lambda * factor1./ (x^2) ', 'x', 'u', 'd', 'factor2', 'lambda', 'factor1');
    end
    
    fprintf('solving %s using newton bisection method... ', se_type); 
     algo_vals(1) = newton_bisect( f, df, dd(1), dd(1) + lambda*norm(v)^2, prox, 1e3 ,  z(1:m), dd(1:m), factor2, lambda, factor1);
    for which_eig = 2:m
         algo_vals(which_eig) = newton_bisect( f, df, dd(which_eig), dd(which_eig - 1), prox, 1e3 ,  z(1:m), dd(1:m), factor2, lambda, factor1);
     end 
    fprintf('Done.\n');

    % calculate eigenvectors
    fprintf('calculating eigenvectors... ');
    diag_org_vals = diag(org_vals);
    for which_eig = 1:m
        temp_diag = 1./(diag_org_vals - algo_vals(which_eig));
        algo_vecs(:, which_eig) = org_vecs * (temp_diag.* z);
    end 
    
    if strcmpi(vf_type, 'none') == 0
                  fprintf('correcting [%s]... ', vf_type);
                 I = find(v~=0);
                 correction_vector = zeros(n,1);
                 for k = 1:n
                     for t=1:size(I,1);
                         if k==I(t)
                             correction_vector(k) = correction_vector(k)  + v(I(t))*(1 - norm(org_vecs(I(t),:))^2);
                         else
                             correction_vector(k) = correction_vector(k)  + v(I(t))*(0 - org_vecs(k,:)*org_vecs(I(t),:)');
                         end
                     end
                 end
                 
                 if strcmpi(vf_type,'TEF') == 1
                     for k = 1:m  
                            algo_vecs(:,k) = algo_vecs(:,k) + correction_vector * (-1/algo_vals(k)); 
                     end
                 elseif strcmpi(vf_type,'CTEF')
                     b = A*correction_vector;
                     for k = 1:m  
                            algo_vecs(:,k) = algo_vecs(:,k) + correction_vector * (-1/algo_vals(k))+  (-1/(algo_vals(k)^2))*b;
                     end
                 end
    end

    
    algo_vecs = normc(algo_vecs);
    
    time = toc;
    
    fprintf('Done.\n');
    
    
end