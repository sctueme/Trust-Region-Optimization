%Pregunta 2.1. 
clear; close all; clc;

[f,grad,Hess] = Camel();
%Un minimo local es (0,0), x0 es cercano a el.
x0 =[0.5;0.5]; 
delta = 1;
fk = f(x0);
g = grad(x0);
B = Hess(x0);
if (norm(x0) <= delta)
    disp('El punto x0 esta dentro de la region de confianza')
end
if (eigs(B,1,'SA') > 0)
    disp('La Hessiana en el punto x0 es definida positiva')
end
pN = -B\g;
pC = pCauchy( B, g, delta );
pDg = pDogLeg( B, g, delta );
M = @(x, y) fk + dot(g,[x;y]) + 0.5*dot([x;y],B*[x;y]);
fsurf(@(x,y) x0(1)+x, @(x,y) x0(2) + y, @(x,y) M(x,y),[-delta,delta,-delta,delta])
colormap gray
hold on
Delta = viscircles(x0', delta);
C = quiver3(x0(1),x0(2),0, pC(1),pC(2),0, 0);
N = quiver3(x0(1),x0(2),0, pN(1),pN(2),0, 0);
Dg = quiver3(x0(1),x0(2),0, pDg(1),pDg(2),0, 0);
legend([Delta, N, C, Dg], { 'Region de confianza', 'Newton','Cauchy', 'Dogleg'});
view(20,20)




