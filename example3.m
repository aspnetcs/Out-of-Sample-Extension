function [ time_ns,  std_ns, time_ms, std_ms] = example3( Ns, m, n, Ms, n_rounds , algo_type)
% Performs runtime comparison between the algorithms in algo_type as function of N amd M
% input:
% (N comparison) Ns = vector of the values of N , i.e., the size of the dataset, m = number of known eigenpairs
% (M comparison) n = size of the dataset, Ms = vector of number of known eigenpairs
% algo_type = cell array, see 'error_and_time.m'
% n_rounds = number of times to perform a single expermient (one N and one
% M
% output:
% time_ns = size(Ns) * (size(algo_type) + 1)  matrix whose columns correspond to mean running time of the algothim 
% in n_rounds experiments and the rows are the corresponding value of N.
% First column is runtime of MATLAB  eigs function
% std_ns = same as time_ns only for the STD of the experiments
% time_ms = size(Ms) * (size(algo_type) + 1) matrix whose columns correspond to mean running time of the algothim 
% in n_rounds experiments and the rows are the corresponding value of M.
% First column is runtime of MATLAB  eigs function
% std_ms = same as time_ms only for the STD of the experiments

    time_ns = zeros(size(Ns,2), size(algo_type,2) +1);
    std_ns = zeros(size(Ns,2), size(algo_type,2) + 1);
    time_ms = zeros(size(Ms,2), size(algo_type,2) +1);
    std_ms = zeros(size(Ms,2), size(algo_type,2) + 1);

    
    % iterations for constant m and various Ns
    for i = 1:size(Ns,2)
        fprintf('n = %d, m = %d\n', Ns(i), m);
        [t, ~, ~] = error_and_time(Ns(i), m, n_rounds, algo_type,  []);
        time_ns(i,:) = mean(t, 1); 
        std_ns(i,:) = std(t, 1); 
    end
                                                                                                                                                                                                                                              
    % iterations for constant n and various Ms
    for i = 1:size(Ms,2)
        fprintf('n = %d, m = %d\n', n, Ms(i));
        [t, ~, ~] = error_and_time(n, Ms(i), n_rounds, algo_type, []);
        time_ms(i,:) = mean(t, 1); 
        std_ms(i,:) = std(t, 1); 
    end

end

