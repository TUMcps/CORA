function out = Gen2_LPV(t,x,u)
global P I

%% State variables -------------------
% Generator 1
delta2   = x(1,:);        %Zustand
omega2   = x(2,:);        %Zustand

%% Uncertain Inputs ----------------------------
V2        = u(2);
theata2   = u(1);

Pe = (I.E2*V2*abs(P.xd2)*cos(angle(P.xd2)+delta2-theata2));

%% Alg. Variables;
dx(1) = P.omegaS*omega2;
dx(2) = (1/P.M2)*(I.Pm2-Pe - (P.D2*omega2));...

out = dx.';