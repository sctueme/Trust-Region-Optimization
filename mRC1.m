function [x, msg] = mRC1( f, x0, itmax )
% Trust region method using the Cauchy point.
%
% In : f ... (handle) function to be optimized
%       x0    ... (vector) initial point
%       itmax ... (natural number) upper bound for number of iterations
%
% Out:  x   ... (vector) last approximation of a stationary point
%       msg ... (string) message that says whether (or not) a minimum was found

    eta = 1e-1;
    tol = 1e-5;
    max_delta = 15e-1;
    delta = 5e-1;
    it = 0;
    xk = x0;
    [~,grad,Hess] = Camel();
    gk = grad(xk);
    Bk = Hess(xk);
    
    while norm(gk, inf) > tol && it < itmax
        pk = pCauchy(Bk, gk, delta);
        coc = (f(xk+pk)-f(xk))/( dot(gk,pk) + 0.5 * dot(pk,Bk*pk));
        if coc < 0.25
            delta = 0.25 * delta;
        elseif coc > 0.75 && (delta - norm(pk)) < eps    
            delta = min(max_delta, 2*delta);
        end
        if coc > eta
            xk = xk + pk;
            gk = grad(xk);
            Bk = Hess(xk);
            it = it + 1;
        end
    end
    x = xk;
    if norm(gk, inf) < tol
        msg = strcat('El metodo converge en ',num2str(it), ' iteraciones');
        [~,p] = chol(Bk);
        if(p == 0)
            msg = strcat(msg, 'El Hessiano es positivo definido, se encontro el minimo local');
        elseif eigs(Bk,1,'SA') >= 0
            msg = strcat(msg, 'Hessiano es positivo semidefinido, se encontro minimo local.');
        end
    else
        msg = strcat('El metodo no converge. Iteraciones: ', num2str(it));
    end
end