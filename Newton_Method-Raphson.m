function Newton_Method()
syms x y 
syms lamda_k

v=[x, y];
global f
f = x^2+y^2+2
%f = sqrt(x^2+(1.2)^2)/2 + sqrt((y-x)^2+(3.7-1.2)^2)/4 + sqrt((4.3-y)^2)/6
g = gradient(f,v)
H = hessian(f,v)
%PASO 1: Escoger un punto inicial X0 y dos parametros de parada E1 y E2
X_k=[5,5]
E1=0.000000001
E2=0.000001
%PASO 2: Calcular el valor del gradiente en X_k
g_Xk=vpa(subs(g,v,X_k))
D=1000
norm(g_Xk)
while (norm(g_Xk)>=E1) & abs(D)>=E1 % PASO 3 
    A=X_k'-g_Xk*lamda_k;
    f_lamda= subs(f,v,A.');
    % Realizar una busqueda unidireccional
    lamda= golden_section(E2,f_lamda,lamda_k)
    H_inv= inv(vpa(subs(H,v,X_k)))  
    if strcmp('FAIL',H_inv)
        break;
    end
    S=H_inv*g_Xk % direction
    vect=lamda*S
    X_k_1=X_k-(vect.');
    f_k=vpa(subs(f,v,X_k));
    f_k_1=vpa(subs(f,v,X_k_1))
    D= (f_k_1-f_k)/f_k;
    X_k=X_k_1;
    g_Xk=vpa(subs(g,v,X_k));
    X_k
    f_k_1
end
end
function [min]= golden_section(epsilon,f,v)

a=0;                            % start of interval
b=10;                            % end of interval
iter= 50;                       % maximum number of iterations
tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
k=0;                            % number of iterations


x1=a+(1-tau)*(b-a);             % computing x values
x2=a+tau*(b-a);

f_x1=vpa(subs(f,v,x1));                     % computing values in x points
f_x2=vpa(subs(f,v,x2));


while ((abs(b-a)>epsilon) && (k<iter))
    k=k+1;
    if(f_x1<f_x2)
        b=x2;
        x2=x1;
        x1=a+(1-tau)*(b-a);
        f_x1=vpa(subs(f,v,x1));                    % computing values in x points
        f_x2=vpa(subs(f,v,x2));
    else
        a=x1;
        x1=x2;
        x2=a+tau*(b-a);
        f_x1=vpa(subs(f,v,x1)) ;                    % computing values in x points
        f_x2=vpa(subs(f,v,x2));
    end
    
    k=k+1;
end


% chooses minimum point
if(f_x1<f_x2)
   min=x1;
   f_x1;
else
   min=x2;
   f_x2;
end
end