function [result]=Lagrange_multiplier()
clear;
syms x y 
syms lamda

f= 4*3.14*(x^2)+2*3.14*y*x
v=[y,x,lamda]
h=2*(3.14*(x^2)*y-4*3.14)

% Plantear la ecuación Lagrangiana L(x,y,lamda)=f(x)+ lamda*sum(h_i)
L=f+ lamda*h
% Resolver el sistema de ecuaciones coyos elemetos son las derivadas
% parciales respecto a cada variable de la funcion lagrangiana igualadas a
% 0
% Resolver el sistema de ecuaciones coyos elemetos son las derivadas
% parciales respecto a cada variable de la funcion lagrangiana igualadas a
% 0
for i=1:length(v)
    system(i)=diff(L,v(i))
end
solutions=[]
solutions = solve(system)
length(solutions.')   
x_values=solutions.x
y_values=solutions.y
lamda_values=solutions.lamda
% hallar el Hessiano ampliado
H = hessian(L,[x,y])
g = gradient(h,[x,y]);
H_a = sym('A', [3 3])
H_a(1,1)=sym('0')
H_a(1,2:length(v))=g
H_a(2:length(v),1)=g.'
H_a(2:length(v),2:length(v))=H
% Evaluar la matriz Hessian Ampliada en las soluciones obtenidas
% Si el determinante es < que cero es mínimo
for i=1:length(lamda_values)
    x_value=x_values(i);
    y_value=y_values(i);
    lamda_value=lamda_values(i);
    H_a_evaluated=subs(H_a, {x, y, lamda}, {x_value,y_value,lamda_value})
    D=det(H_a_evaluated)
    if D<0
        minima=[x_value,y_value,lamda_value]
    else
      maxima=[x_value,y_value,lamda_value]
    end    
end
end