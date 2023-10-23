function f = dynamics_quad(x,u)
% dynamics_quad - neural-network controlled quadrotor [1]
%
% Syntax:
%    f = dynamics_quad(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
% 
% Outputs:
%    f - time-derivate of the system state
%
% Reference:
%   [1] Beard, R.: Quadrotor dynamics and control rev 0.1 (2008)

% Authors:       Tobias Ladner
% Written:       20-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parameter
g = 9.81;
m = 1.4;
Jx = 0.054;
Jy = 0.054;
Jz = 0.104;
tau = 0;

% extract variables for better readability
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);
x5 = x(5);
x6 = x(6);
x7 = x(7);
x8 = x(8);
x9 = x(9);
x10 = x(10);
x11 = x(11);
x12 = x(12);

u1 = u(1);
u2 = u(2);
u3 = u(3);

% dynamic equations


f = [
     cos(x8)*cos(x9)*x4 + (sin(x7)*sin(x8)*cos(x9) - cos(x7)*sin(x9))*x5 ...
        + (cos(x7)*sin(x8)*cos(x9) + sin(x7)*sin(x9))*x6;
     cos(x8)*sin(x9)*x4 + (sin(x7)*sin(x8)*sin(x9) + cos(x7)*cos(x9))*x5 ...
        + (cos(x7)*sin(x8)*sin(x9) - sin(x7)*cos(x9))*x6;
     sin(x8)*x4 - sin(x7)*cos(x8)*x5 - cos(x7)*cos(x8)*x6;
     x12*x5 - x11*x6 - g*sin(x8);
     x10*x6 - x12*x4 + g*cos(x8)*sin(x7);
     x11*x4 - x10*x5 + g*cos(x8)*cos(x7) - g - u1/m;
     x10 + sin(x7)*tan(x8)*x11 + cos(x7)*tan(x8)*x12;
     cos(x7)*x11 + sin(x7)*x12;
     sin(x7)/cos(x8) * x11 - sin(x7)*x12;
     (Jy-Jz)/Jx * x11*x12 + 1/Jx * u2;
     (Jz-Jx)/Jy * x10*x12 + 1/Jy * u3;
     (Jx-Jy)/Jz * x10*x11 + 1/Jz * tau
];
     
end

% ------------------------------ END OF CODE ------------------------------
