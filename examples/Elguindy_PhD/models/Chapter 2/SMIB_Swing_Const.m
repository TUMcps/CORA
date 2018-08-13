function out = SMIB_Swing_Const(t,x,y,u)
global P I

if strcmp(P.mode,'normal')
    xs = P.xs;
elseif strcmp(P.mode,'fault')
    xs = P.xs3;
end


%% State variables -------------------
delta1   = x(1);        %Zustand


%% Alg. varibles ---------------------
P1      = y(1,:);        % active power at bus 1 (alg. variable)
Q1      = y(2,:);        % reactive power at bus 1 (alg. variable)
Theta1  = y(3,:);       % phase angle at bus 1 (alg. variable)
V1      = y(4,:);       % voltage at bus 1 (alg. variable)

Pe1 = (I.E1*V1*abs(P.xd1)*cos(angle(P.xd1)+delta1-Theta1));

%% Alg. Variables;
dy(1) = P1 - Pe1;
dy(2) = Q1 + (I.E1*V1*abs(P.xd1)*sin(angle(P.xd1)+delta1-Theta1)) - (V1^2*abs(P.xd1)*sin(angle(P.xd1)));
dy(3) = V1*P.v2/xs*sin(Theta1-P.theta2)-P1;
dy(4) = (V1*V1/xs)-V1*P.v2/xs*cos(Theta1-P.theta2)-Q1;


out = dy.';