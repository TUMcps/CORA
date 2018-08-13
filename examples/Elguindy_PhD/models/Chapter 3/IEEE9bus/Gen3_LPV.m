function out = Gen3_LPV(t,x,u)
global P I

%% State variables -------------------
% Generator 1
delta3   = x(1,:);        %Zustand
omega3   = x(2,:);        %Zustand

%% Uncertain Inputs ----------------------------
V3        = u(2);
theata3   = u(1);

Pe = (I.E3*V3*abs(P.xd3)*cos(angle(P.xd3)+delta3-theata3));

%% Alg. Variables;
dx(1) = P.omegaS*omega3;
dx(2) = (1/P.M3)*(I.Pm3-Pe - (P.D3*omega3));...

out = dx.';