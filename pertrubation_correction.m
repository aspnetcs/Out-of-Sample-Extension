function [algo_vecs_cor, algo_vals_cor] = pertrubation_correction(algo_vecs, algo_vals, C)

% implements the perturbation correction algorithm with correction matrix C
        
        m = size(algo_vecs, 2);
        algo_vals_cor = zeros(m,1);
        algo_vecs_cor = zeros(size(algo_vecs));
        
        algo_vals_cor = diag(diag(algo_vals) + algo_vecs' * C * algo_vecs);


        
        for i = 1:m
            algo_vecs_cor(:,i) = algo_vecs(:,i);
            mul = (C*algo_vecs(:,i));
            for j = 1:m
                if j~=i
                    algo_vecs_cor(:,i)  = algo_vecs_cor(:,i)  + ((algo_vecs(:,j)'* mul)/(algo_vals(i) - algo_vals(j))) * algo_vecs(:,j);
                end
                
            end
        end
      
        algo_vecs_cor = normc(algo_vecs_cor);

end

