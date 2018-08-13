function [dx]=tank6paramEq(t,x,u,p)

k=0.015;

k2=0.01;

g=9.81; 


%dx(1,1)=u(1)+0.1+k2*exp(4-x(6))-p(1)*sqrt(2*g)*sqrt(x(1)); %tank 1
dx(1,1)=u(1)+0.1+k2*(4-x(6))-p(1)*sqrt(2*g)*sqrt(x(1)); %tank 1
dx(2,1)=p(1)*sqrt(2*g)*(sqrt(x(1))-sqrt(x(2))); %tank 2
dx(3,1)=p(1)*sqrt(2*g)*(sqrt(x(2))-sqrt(x(3))); %tank 3
dx(4,1)=p(1)*sqrt(2*g)*(sqrt(x(3))-sqrt(x(4))); %tank 4
dx(5,1)=p(1)*sqrt(2*g)*(sqrt(x(4))-sqrt(x(5))); %tank 5
dx(6,1)=p(1)*sqrt(2*g)*(sqrt(x(5))-sqrt(x(6))); %tank 6