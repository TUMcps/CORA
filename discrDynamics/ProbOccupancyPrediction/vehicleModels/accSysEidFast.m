function dx =accSysEidFast(x,u)

aMax=7;
c1=7.32;

% dx(1,1)=1e-3*x(1)+x(2); %position
% dx(2,1)=1e-3*x(2)+aMax*c1/x(2)*u; %speed

dx(1,1)=x(2); %position
dx(2,1)=aMax*c1/x(2)*u; %speed