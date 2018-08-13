function out = Gen1_30(t,x,u)
global P I

%% State variables -------------------
% Generator 1
delta1   = x(1,:);        %Zustand
omega1   = x(2,:);        %Zustand

%% Uncertain Inputs ----------------------------
V1        = u(1);
theata1   = I.theta1;

Pe = (I.E1*V1*abs(P.xd1)*cos(angle(P.xd1)+delta1-theata1));

%% Alg. Variables;
dx(1) = P.omegaS*omega1;
dx(2) = (1/P.M1)*(I.Pm1-Pe - (P.D1*omega1));...

out = dx.';