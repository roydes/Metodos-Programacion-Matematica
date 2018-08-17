%**************** HOOKE-JEEVES METHOD ***********************
function [x]= Hooke_Jeeves_Method()
%--------------STEP-1------------------------
xk=[0,0];
delta=[0.5,0.5]
alpha=2;
epsilon= 0.001
x=[0,0]
while true
  %--------------STEP-2------------------------
  % Exploratory Move
   x= Exploratory_Move(xk,delta)
  %--------------STEP-3------------------------
    if x~=xk
        xk1=x  
        while ObjectiveFunction(xk1)< ObjectiveFunction(xk)
  %--------------STEP-4------------------------
  %Pattern Move
            xp=xk1+(xk1-xk)
 %--------------STEP-5------------------------
  % Exploratory Move
            x=Exploratory_Move(xp,delta)
            xk=xk1
            xk1=x
        end
    else if norm(delta)< epsilon
            break;
        else delta=delta/alpha
        end
        
    end
    
end
end


function [result]= Exploratory_Move(xk,delta)
n=length(xk)
x=xk;
for i = 1:n
   xtemp2=x
   xtemp3=x
   xtemp2(i)=xtemp2(i)+delta(i);
   xtemp3(i)=xtemp3(i)-delta(i);
   
    f1= ObjectiveFunction(x)
    f2= ObjectiveFunction(xtemp2)
    f3= ObjectiveFunction(xtemp3)
    
 
    fmin=min(f1,min(f2,f3))
    if fmin==f1
        x=x
    elseif fmin==f2
        x=xtemp2
    else x=xtemp3
    end
    
end
result=x

end






function y =  ObjectiveFunction(x)
%Override function
y = (x(1)^2+x(2)-11)^2+(x(1)+x(2)^2-7)^2;
%y=(1-x(1))^2+100*(x(2)-x(1)^2)^2;
end