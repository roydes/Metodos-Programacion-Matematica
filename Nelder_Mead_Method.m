
%**************** NELDER-MEAD METHOD ***********************
function [minimum]= Nelder_Mead_Method()
%STEP-1
%Choosing gamma>1, 0<beta<1
gamma=1.5;
%input('Input a value for gamma')
beta=0.5;
%input('Input a value for beta')
epsilon=0.00001;
%input('Input a value for epsilon')
%Creating an initial simplex
n= 3;
x=rand(n,n-1);
Q=10000;
minimum=10000;
%STEP-2
while Q > epsilon
    %Sorting points by the objective function value
    for i = 1:n
        for j = 1:n-i
            if (ObjectiveFunction(x(j,:)) < ObjectiveFunction(x(j+1,:)))
                p= x(j,:);
                x(j,:) = x(j+1,:);
                x(j+1,:) = p;
            end
        end
    end
    
    %Highest function evaluation-Worst point
    xh=x(1,:)
    %Next to the worst point
    xg=x(2,:)
    %Lowest-function evaluation-best point
    xl=x(n,:)
    %Calculating the centroid
    xc=sum(x(2:n,1:end))/(n-1);
    %STEP 3
    %Reflection
    xr=2*xc-xh;
    xnew=xr;
    minimum=xl
    if ObjectiveFunction(xr)<ObjectiveFunction(xl)
        %Expansion
        xnew= (1+gamma)*xc-gamma*xh;
    elseif ObjectiveFunction(xr)>= ObjectiveFunction(xh)
        %Contraction
        xnew= (1-beta)*xc+beta*xh;
    elseif ObjectiveFunction(xg)<ObjectiveFunction(xr)&& ObjectiveFunction(xr)<ObjectiveFunction(xh)
        %Contraction
        xnew= (1+beta)*xc-beta*xh;
    end
    s=0;
    x(1,:)=xnew
    
    for i = 1:n
        s=s+(ObjectiveFunction(x(i,:))-ObjectiveFunction(xc))^2;
    end
    Q=sqrt(s/n)
    
end

end

function y =  ObjectiveFunction(x)
%Override function
%y = (x(1)^2+x(2)-11)^2+(x(1)+x(2)^2-7)^2;
%y=(1-x(1))^2+100*(x(2)-x(1)^2)^2;
y= sin(x(1)+x(2))+(x(1)-x(2))^2-1.5*x(1)+2.5*x(2)+1;
end



