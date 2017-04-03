function [ L0, L1, min_k] = build_graph_laplacian( X, IDX, epsilon )
%calculates the graph Laplcian matrix of X, with neighbours as indicated in
%IDX and paramter epsilon for the Gaussian kernel. L0 is the graph Laplcian
%with the first node disconnected, and L1 is the graph Laplican after it
%was connected.

        n_samples = size(X,1);
        k = size(IDX,2);
        fprintf('calculating graph Laplacian matrices... '); 
        W1 = zeros(n_samples,n_samples);
        for m = 1:n_samples
            for n = 1:k
                 curr = IDX(m,n);
                 xm = X(m,:); xn = X(n,:);
                 value = exp(-norm(xm - xn)^2/epsilon);
                 W1(m,curr) = value;
                 W1(curr,m) = value;
            end
        end
        W1 = sparse(W1);
        W0 = W1; W0(1,:) = 0; W0(:,1) = 0; W0(1,1) = 1;
        Dinv0 = sparse(diag(1./sum(W0,2))); L0  = (Dinv0.^(0.5))*(W0)* (Dinv0.^(0.5));
        clear W0;clear Dinv0;
        D1 = sparse(diag(sum(W1,2))); Dinv1 = sparse(diag(1./sum(W1,2))); L1  = (Dinv1.^(0.5))*(W1)* (Dinv1.^(0.5));
        min_k= min(diag(D1(2:end, 2:end))); 
        clear W1; clear D1; clear Dinv1;
        
        fprintf('done.\n');  
end

