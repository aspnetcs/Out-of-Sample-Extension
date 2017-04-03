function [L0, L1] = huge_example(N, k)

    % DO L0
    I = zeros((N-1) * (k +1) + 1, 1);
    J = zeros((N-1) * (k +1) + 1, 1);
    idx = 1;
    for i=1:N
        if (i == 1)
                I(idx) = 1;
                J(idx) = 1;
                idx = idx + 1;
        else
            for j = (i - k/2):(i + k/2)
                         I(idx) = i;

                         if (j < 1)
                             place = j + N;
                         elseif (j > N)
                             place = j - N;
                         elseif( j== 1)  % first vertex, ignore
                             place = 1;
                         else
                            place = j;
                         end
                         if (place~=1)
                                J(idx) = place;
                                idx = idx + 1;
                         end
            end
        end
    end
    I(I == 0) = []; J(J == 0) = [];
    L0 = sparse(I,J, 1/(k+1)); L0(1,1) = 1;

    % DO L1
    I = zeros((N-1) * (k +1) + 1, 1);
    J = zeros((N-1) * (k +1) + 1, 1);
    idx = 1;
    for i=1:N
            for j = (i - k/2):(i + k/2)
                         I(idx) = i;

                         if (j < 1)
                             place = j + N;
                         elseif (j > N)
                             place = j - N;
                         elseif( j== 1)  % first vertex, ignore
                             place = 1;
                         else
                            place = j;
                         end
                         J(idx) = place;
                          idx = idx + 1;
            end
    end
    I(I == 0) = []; J(J == 0) = [];
    L1 = sparse(I,J, 1/(k+1)); 
end