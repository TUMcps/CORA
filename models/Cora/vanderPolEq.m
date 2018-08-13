function [dx]=vanderPolEq(t,x,u)

mu=1;
%mu=0.3;

dx(1,1)=x(2);
dx(2,1)=mu*(1-x(1)^2)*x(2)-x(1)+u(1);
