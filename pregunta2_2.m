%Pregunta 2.2
clear; close all; clc;
[f,grad,hess] = Camel();
%Un minimo local es (0,0), Elegimos x0 y vemos que sea lejano a (0,0)
x0 =[1;3];
radioMax = 1.5;
g = grad(x0);
B = hess(x0);
if (norm(x0) > radioMax)
    disp('El punto x0 esta fuera de la region de confianza')
end
if (eigs(B,1,'SA') <= 0)
    disp('La Hessiana en el punto x0 no es definida positiva')
end

for itMax = [100:100:1000]
    disp(strcat('Iteraciones: ',num2str(itMax)))
    tol = 10^(-5);
    x_opt = [0;0];
    radioMax = 1.5;
    disp('Con MRC1')
    [x,msg] = mRC1(f, x0, itMax)
    error = norm(x-x_opt)
    disp('Con MRC2')
    [x,msg] = mRC2(f, x0, itMax)
    error = norm(x-x_opt)
end

