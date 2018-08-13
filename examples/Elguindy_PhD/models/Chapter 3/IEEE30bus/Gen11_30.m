function out = Gen11_30(t,x,u)
global P I

%% State variables -------------------
% Generator 1
delta11   = x(1,:);        %Zustand
omega11   = x(2,:);        %Zustand

%% Uncertain Inputs ----------------------------
V11        = u(2);
theata11   = u(1);

Pe = (I.E11*V11*abs(P.xd11)*cos(angle(P.xd11)+delta11-theata11));

%% Alg. Variables;
dx(1) = P.omegaS*omega11;
dx(2) = (1/P.M11)*(I.Pm11-Pe - (P.D11*omega11));...

out = dx.';