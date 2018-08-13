function dx = Model_2D(t,x,u)

x1 = x(1);
x2 = x(2);

dx(1,1) = -2*x1 + x2 + x1^3 + x2^5;
dx(2,1)= -x1-x2+x1^2*x2^3;  