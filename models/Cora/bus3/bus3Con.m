function f = bus3Con(x,y,u)
% bus2Con - constraint function for a 3-bus power system
%
% Syntax:  
%    f = bus3Con(x,y,u)
%
% Inputs:
%    x - state vector
%    y - vector of algebraic variables
%    u - input vector
%
% Outputs:
%    f - vector storing the constraint viaolation

% Author:       Niklas Kochdumper
% Written:      19-June-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%x(1) = omega
%x(2) = T_m

%y(1) = E
%y(2) = V_2
%y(3) = V_3
%y(4) = Theta_1 - delta
%y(5) = Theta_2 - delta
%y(6) = Theta_3 - delta

%u(1) = P_c
%u(2) = P_w

%slack bus voltage
V_1 = 1;


%admittances
Z_m = 0.2j;
Z_13 = 0.1j;
Z_23 = 0.15j;

Y_m = abs(1/Z_m);
Y_13 = abs(1/Z_13);
Y_23 = abs(1/Z_23);

%parameters
P_3 = 1;
Q_3 = 0.5;
Q_w = 0; %reactive wind power


%constraints
f(1,1) = Y_m*y(1)*V_1*sin(-y(4)) - Y_13*V_1*y(3)*sin(y(4) - y(6)); %eq (21)
f(2,1) = u(2) - Y_23*y(2)*y(3)*sin(y(5) - y(6)); %eq (23)
f(3,1) = -P_3 - Y_13*V_1*y(3)*sin(y(6) - y(4)) - Y_23*y(2)*y(3)*sin(y(6) - y(5)); %eq (25)
f(4,1) = Y_m*y(1)*V_1*cos(-y(4)) - Y_m*V_1^2 + Y_13*V_1*y(3)*cos(y(4) - y(6)) - (Y_13 + Y_m)*V_1^2; %eq (22)
f(5,1) = Q_w + Y_23*y(2)*y(3)*cos(y(5) - y(6)) - Y_23*y(2)^2; %eq (24)
f(6,1) = -Q_3 + Y_13*V_1*y(3)*cos(y(6) - y(4)) + Y_23*y(2)*y(3)*cos(y(6) - y(5)) - (Y_13 + Y_23)*y(3)^2; %eq (26)

%------------- END OF CODE --------------