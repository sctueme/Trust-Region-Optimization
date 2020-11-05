%Funcion Three-Hump Camel
function [f,Grad,Hess] = Camel()

f = @(x) 2*(x(1))^2 -1.05*x(1)^4 + (x(1)^6)/6 + x(1)*x(2) + x(2)^2;

Grad = @(x) [(x(1)^5) - 4.2*x(1)^3 + 4*x(1) + x(2);x(1) + 2*x(2)];

Hess = @(x) [4 - 12.6*x(1)^2 + 5*x(1)^4, 1;1, 2];

end

