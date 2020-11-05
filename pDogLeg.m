function [p] = pDogLeg( B, g, delta )
    % In : B ... an s.p.d. matrix that approximates the hessian of f in xk
    %       g     ... (vector) gradient of  f  in  xk
    %       delta ... trust region radius
    %
    % Out:  p     ... The  dogleg  point

    pC = pCauchy(B,g,delta);
    pN = -B\g; % <=> pN = -inv(B)*g 
    if(norm(pN) <= delta)
       p = pN;
    else
        if(norm(pC) > delta)
            p = -(delta/norm(g))*g;
        else
            alphas = roots([(norm(pN-pC))^2,2*dot((pN-pC),pC), (norm(pC))^2 - (delta^2)]);
            alpha = max(alphas);
            p = pC+ alpha*(pN-pC);
        end

    end
end