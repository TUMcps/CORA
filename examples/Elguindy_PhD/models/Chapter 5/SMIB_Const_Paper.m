function out = SMIB_Const_Paper(t,x,y,u)
global P

if strcmp(P.mode,'normal')
    xs = P.xs;
elseif strcmp(P.mode,'fault')
    xs = P.xs3;
end


%% State variables -------------------
delta   = x(1);        %Zustand
omega   = x(2);        %Zustand
Eq_     = x(3);        %Zustand 

%% Inputs ----------------------------
deltarl  = u(1);
omegarl  = u(2);
Eq_rl    = u(3);

%% Alg. varibles ---------------------
ed      = y(1);        %alg. Variable
eq      = y(2);        %alg. Variable
id      = y(3);        %alg. Variable
iq      = y(4);        %alg. Variable
P1      = y(5);        %alg. Variable
Q1      = y(6);        %alg. Variable
Theta1  = y(7);        %alg. Variable
V1      = y(8);        %alg. Variable

%% Alg. Variables;
dy(1) = -ed+V1*sin(delta-Theta1);
dy(2) = -eq+V1*cos(delta-Theta1);
dy(3) = -ed+P.xq*iq-P.ra*id;
dy(4) = -eq+(Eq_)-P.ra*iq-P.x_d*id;
dy(5) = -P1+ed*id+eq*iq;
dy(6) = -Q1+eq*id-ed*iq;
dy(7) = V1*P.v2/xs*sin(Theta1-P.theta2)-P1;
dy(8) = V1*V1/xs-V1*P.v2/xs*cos(Theta1-P.theta2)-Q1;

out = dy.';