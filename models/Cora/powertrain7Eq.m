function f=powertrain7Eq(t,x,u,p)

%control
v = p.k_K*(p.i*x(4) - x(7)) ...
    + p.k_KD*(p.i*u(1) - 1/p.J_m*(x(2) - 1/p.i*p.k*(x(1) - p.alpha) - p.b_m*x(7))) ...
    + p.k_KI*(p.i*x(3) - p.i*(x(1) + x(5))) ...
    + 0*1/p.i*p.J_l*u(1);

%plant model
f(1,1) = 1/p.i*x(7) - x(6); %Theta_d
f(2,1) = (v - x(2))/p.tau_eng; %T_m
f(3,1) = x(4); %Theta_ref
f(4,1) = u(1); %\dot{Theta}_ref
f(5,1) = x(6); %Theta_l
f(6,1) = 1/p.J_l*(p.k*(x(1) - p.alpha) - u(2) - p.b_l*x(6)); %\dot{Theta}_l
f(7,1) = 1/p.J_m*(x(2) - 1/p.i*p.k*(x(1) - p.alpha) - p.b_m*x(7)); %\dot{Theta}_m
