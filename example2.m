function [errors] = example2(n, m, mu_sizes)
% returns the errors of a synthetic rank one update by several algorithms
% and plot log-error vs log-mu for TSE,CTSE, TCVF1 and TCFV2
% input:
% n - size of the matrix
% m - number of known eigenpairs
% mu_size - vector contating the size of the unkown eigenvalues
% output:
% errors: size(mu_sizes) * 6 matrix, whose columns  are: mu values, TSE error, CTSE error, TVF error, TCVF1 error and TCFV2 error 

    A = randn(n); A = A'*A; A = A /norm(A); [Q,D] = eig(A); [~,P] = sort(diag(D), 'descend'); Q = Q(:,P); D = D(P,P);
    
    v = 1 + randn(n,1); v = v/norm(v);
    errors_tse = zeros(size(mu_sizes,2), 1);
    errors_ctse = zeros(size(mu_sizes,2), 1);
    errors_tse_vecs = zeros(size(mu_sizes,2), 1);
    errors_ctse_vecs = zeros(size(mu_sizes,2), 1);
    errors_ctse_vecs_corr1 = zeros(size(mu_sizes,2), 1);
    errors_ctse_vecs_corr0 = zeros(size(mu_sizes,2), 1);
    
    for i = 1:size(mu_sizes,2)
        
        for j=(m+1):n
                D(j,j) = mu_sizes(i);
        end
        [~,P] = sort(diag(D), 'descend'); Q = Q(:,P); D = D(P,P);
        A = Q*D*Q';

        [real_vecs, real_vals] = eigs(A + (v*v'), m);
        [~,P] = sort(diag(real_vals), 'descend'); real_vecs = real_vecs(:,P); real_vals = real_vals(P,P);
        [old_vecs, old_vals] = eigs(A, m);
        [~,P] = sort(diag(old_vals), 'descend'); old_vecs = old_vecs(:,P); old_vals = old_vals(P,P);
        
        [tse_eigvecs, tse_eigvals, ~] = update_eigenspectrum(A, 1, v, old_vecs, old_vals, 0, '11');
        [ctse_eigvecs, ctse_eigvals, ~] = update_eigenspectrum(A, 1, v, old_vecs, old_vals, 0, '22');

        [ctse_vecs_corr1 , ~, ~] = update_eigenspectrum(A, 1, v, old_vecs, old_vals, 0, '21');
        [ctse_vecs_corr0 , ~, ~] = update_eigenspectrum(A, 1, v, old_vecs, old_vals, 0, '20');
        
        errors_tse(i) = mean(abs(tse_eigvals - diag(real_vals)));
        errors_ctse(i) = mean(abs(ctse_eigvals - diag(real_vals)));
        errors_tse_vecs(i) = mean((min(sum((real_vecs - tse_eigvecs).^2,1), sum((real_vecs+ tse_eigvecs).^2,1))).^0.5);
        errors_ctse_vecs(i) = mean((min(sum((real_vecs - ctse_eigvecs).^2,1), sum((real_vecs+ ctse_eigvecs).^2,1))).^0.5);
        errors_ctse_vecs_corr1(i) = mean((min(sum((real_vecs - ctse_vecs_corr1).^2,1), sum((real_vecs+ ctse_vecs_corr1).^2,1))).^0.5);
        errors_ctse_vecs_corr0(i) = mean((min(sum((real_vecs - ctse_vecs_corr0).^2,1), sum((real_vecs+ ctse_vecs_corr0).^2,1))).^0.5);
    end

    figure;
    hold on
    scatter(log2(mu_sizes), log2(errors_tse), 'r', '+');
    scatter(log2(mu_sizes), log2(errors_ctse), 'b', 'o');
    xlabel('log2(mu)');
    ylabel('log2(mean error)');
    title('eigenvalues mean error');
    legend('TSE', 'CTSE');
    hold off
    
    figure;
    hold on
    scatter(log2(mu_sizes), log2(errors_ctse_vecs_corr1), 'r', '+');
    scatter(log2(mu_sizes), log2(errors_ctse_vecs), 'b', 'o');
    xlabel('log2(mu)');
    ylabel('log2(mean error)');
    title('eigenvectors mean error');
    legend('correction 1', 'correction 2');
    hold off
    
    errors = [mu_sizes', errors_tse, errors_ctse, errors_ctse_vecs_corr0, errors_ctse_vecs_corr1, errors_ctse_vecs];
    
end

