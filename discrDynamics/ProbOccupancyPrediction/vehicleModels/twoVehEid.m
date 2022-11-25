function [dx]=twoVehEid(t,x,u)

%x1: position following vehicle
%x2: position leading vehicle
%x3: velocity following vehicle
%x4: velocity leading vehicle


aMax=7;
c1=7.32;

%position updates
dx(1,1)=x(3); %following vehicle
dx(2,1)=x(4); %leading vehicle

%velocity update for following vehicle
if (x(3)<=0)
    dx(3,1)=0;
else
    if (x(3)<c1) | (u(1)<0)
        dx(3,1)=aMax*u(1); 
    else
        dx(3,1)=aMax*c1/x(3)*u(1); 
    end
end

%velocity update for leading vehicle
if (x(4)<=0)
    dx(4,1)=0;
else
    if (x(4)<c1) | (u(2)<0)
        dx(4,1)=aMax*u(2); 
    else
        dx(4,1)=aMax*c1/x(4)*u(2); 
    end
end