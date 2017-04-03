function [ x] = newton_bisect( f, df, a, b, tol, nmax ,  u, d, factor2, lambda, factor1)
% implements newton/bisection method for solving equation f

    if abs(f(a,u, d, factor2, lambda, factor1)) < tol
        x = a;
        return;
    elseif abs(f(b,u, d, factor2, lambda, factor1)) < tol
        x = b;
        return;
    end
    
    a = a + 10*eps;
    b = b - 10*eps;
    
    x(1) = 0.5*(a + b);
    ex(1) = abs(x(1)- a);
    
    k = 2;
    
    while (k <= nmax)
        fdx = (f(x(k-1),u, d, factor2, lambda, factor1)/df(x(k-1),u, d, factor2, lambda, factor1));
        if (abs(b - a) < tol) || (abs(fdx) < tol)
        %if abs(f(x(k-1),u, d, A, lambda)) < tol
            x = x(end);
            return;
        else
            
            x(k) = x(k-1) - fdx;
            
            if ((x(k) < a) || (x(k) > b)) % out of [a,b], do bisection
                
                x(k) = 0.5 * (a+b);
                
                if (f(a, u,d, factor2, lambda, factor1) * f(x(k),u,d, factor2, lambda, factor1) < 0)
                    b = x(k);
                else
                    a = x(k);
                end   
            end
            ex(k) = abs(x(k)-x(k-1));
        end
        k = k + 1;
    end

    
    x = x(end);

end