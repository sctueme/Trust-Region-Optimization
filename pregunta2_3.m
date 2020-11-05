%Pregunta 2.3
clear;   close all;   clc;
L = @(x,y) 2*(x.^2) -1.05*(x.^4) + (x.^6/6) + (x.*y) + (y.^2);
stepsize =  0.01; 
[X,Y] = meshgrid(-2:stepsize:2);
z = L(X,Y);
niveles = [0.1, -5:5];
contour(X,Y,z, niveles)

axis equal

hold on
x0 = [1;3];
x = [2:-0.2:0];
iterRC1 = 200;
iterRC2 = 10;

Cauchy = zeros(iterRC1, 2);
Dogleg = zeros(iterRC2, 2);
[f,grad,hess] = Camel();
xk = x0;
for i = 1:iterRC1
    Cauchy(i,:) = xk;
    [xk, ~] = mRC1(f, xk, 1);
end

xk = x0;
for i = 1:iterRC2
    Dogleg(i,:) = xk;
    [xk, ~] = mRC2(f, xk, 1);
end

trayC = plot(Cauchy(:, 1), Cauchy(:, 2),'--d');
trayD = plot(Dogleg(:, 1), Dogleg(:, 2),'--+');
legend([trayC, trayD], { 'Cauchy', 'Dogleg'});
view(20,20)
hold off
