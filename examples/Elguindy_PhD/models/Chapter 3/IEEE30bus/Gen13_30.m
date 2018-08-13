function out = Gen13_30(t,x,u)
global P I

%% State variables -------------------
% Generator 1
delta13   = x(1,:);        %Zustand
omega13   = x(2,:);        %Zustand

%% Uncertain Inputs ----------------------------
V13        = u(2);
theata13   = u(1);

Pe = (I.E13*V13*abs(P.xd13)*cos(angle(P.xd13)+delta13-theata13));

%% Alg. Variables;
dx(1) = P.omegaS*omega13;
dx(2) = (1/P.M13)*(I.Pm13-Pe - (P.D13*omega13));...

out = dx.';