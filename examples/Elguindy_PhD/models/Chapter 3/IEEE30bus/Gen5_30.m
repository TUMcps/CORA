function out = Gen5_30(t,x,u)
global P I

%% State variables -------------------
% Generator 1
delta5   = x(1,:);        %Zustand
omega5   = x(2,:);        %Zustand

%% Uncertain Inputs ----------------------------
V5        = u(2);
theata5   = u(1);

Pe = (I.E5*V5*abs(P.xd5)*cos(angle(P.xd5)+delta5-theata5));

%% Alg. Variables;
dx(1) = P.omegaS*omega5;
dx(2) = (1/P.M5)*(I.Pm5-Pe - (P.D5*omega5));...

out = dx.';