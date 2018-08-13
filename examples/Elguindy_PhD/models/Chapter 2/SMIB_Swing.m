function out = SMIB_Swing(t,x,y,u)
global P I

%% State variables -------------------
delta1   = x(1);        %Zustand
omega1   = x(2);        %Zustand


%% Alg. varibles ---------------------
P1      = y(1,:);        % active power at bus 1 (alg. variable)
Q1      = y(2,:);        % reactive power at bus 1 (alg. variable)
Theta1  = y(3,:);       % phase angle at bus 1 (alg. variable)
V1      = y(4,:);       % voltage at bus 1 (alg. variable)

Pe1 = (I.E1*V1*abs(P.xd1)*cos(angle(P.xd1)+delta1-Theta1));


%% State Varaibles
dx(1) = P.omegaS*omega1;
dx(2) = (1/P.M1)*(I.Pm1-Pe1 - (P.D1*omega1));

out = dx.';