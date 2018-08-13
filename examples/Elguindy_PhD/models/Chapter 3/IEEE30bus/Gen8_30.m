function out = Gen8_30(t,x,u)
global P I

%% State variables -------------------
% Generator 1
delta8   = x(1,:);        %Zustand
omega8   = x(2,:);        %Zustand

%% Uncertain Inputs ----------------------------
V8        = u(2);
theata8   = u(1);

Pe = (I.E8*V8*abs(P.xd8)*cos(angle(P.xd8)+delta8-theata8));

%% Alg. Variables;
dx(1) = P.omegaS*omega8;
dx(2) = (1/P.M8)*(I.Pm8-Pe - (P.D8*omega8));...

out = dx.';