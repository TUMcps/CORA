function [dx]=fiveDimSysEq(t,x,u)

dx(1,1) = -x(1) - 4*x(2) + u(1); %dimension 1
dx(2,1) = 4*x(1) - x(2) + u(2); %dimension 2
dx(3,1) = -3*x(3) + x(4) + u(3); %dimension 3
dx(4,1) = -x(3) -3*x(4) + u(4); %dimension 4
dx(5,1) = -2*x(5) + u(5); %dimension 5