function dx = quadrocopterControlledEq(x,u)
% quadrocopterControlledEq - dynamic equations for the quadrocopter 
%                            benchmark (see Eq (16) - (19) in [1])
%
% Syntax:  
%    dx = quadrocopterControlledEq(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
%
% Outputs:
%    dx - time-derivate of the system state
% 
% References:
%   [1] R. Beard, "Quadrotor Dynamics and Control", Tech Report Bringham
%       Young University

% Author:        Ahmed El-Guindy, Matthias Althoff, Mark Wetzlinger
% Written:       22-May-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% x_1 = p_n
% x_2 = p_e
% x_3 = h

% x_4 = u
% x_5 = v
% x_6 = w

% x_7 = phi
% x_8 = theta
% x_9 = psi

% x_10 = p
% x_11 = q
% x_12 = r

% u_1 desired height

% paremeters
g = 9.81; %[m/s^2], gravity constant
R = 0.1;
l = 0.5;
M_rotor = 0.1;
M = 1;
m = M + 4*M_rotor;

% auxiliary parameters
J_x = 2*M*R^2/5 + 2*l^2*M_rotor;
J_y = J_x;
J_z = 2*M*R^2/5 + 4*l^2*M_rotor;

% height control (PD)
F = m*g - 10*(x(3) - u(1)) + 3*x(6);

% desired 

% roll control (PD)
tau_phi = -(x(7) - u(2)) - x(10);

% pitch control (PD)
tau_theta = -(x(8) - u(3)) - x(11);

% heading is uncontrolled:
tau_psi = 0;

dx(1,1) = cos(x(8))*cos(x(9))*x(4) + (sin(x(7))*sin(x(8))*cos(x(9)) - cos(x(7))*sin(x(9)))*x(5) + (cos(x(7))*sin(x(8))*cos(x(9)) + sin(x(7))*sin(x(9)))*x(6);
dx(2,1) = cos(x(8))*sin(x(9))*x(4) + (sin(x(7))*sin(x(8))*sin(x(9)) + cos(x(7))*cos(x(9)))*x(5) + (cos(x(7))*sin(x(8))*sin(x(9)) - sin(x(7))*cos(x(9)))*x(6);
dx(3,1) = sin(x(8))*x(4) - sin(x(7))*cos(x(8))*x(5) - cos(x(7))*cos(x(8))*x(6);

dx(4,1) = x(12)*x(5) - x(11)*x(6) - g*sin(x(8));
dx(5,1) = x(10)*x(6) - x(12)*x(4) + g*cos(x(8))*sin(x(7));
dx(6,1) = x(11)*x(4) - x(10)*x(5) + g*cos(x(8))*cos(x(7)) - F/m;

dx(7,1) = x(10) + (sin(x(7))*tan(x(8)))*x(11) + (cos(x(7))*tan(x(8)))*x(12);
dx(8,1) = cos(x(7))*x(11) - sin(x(7))*x(12);
dx(9,1) = sin(x(7))/cos(x(8))*x(11) + cos(x(7))/cos(x(8))*x(12);

dx(10,1) = (J_y - J_z)/J_x*x(11)*x(12) + 1/J_x*tau_phi;
dx(11,1) = (J_z - J_x)/J_y*x(10)*x(12) + 1/J_y*tau_theta;
dx(12,1) = (J_x - J_y)/J_z*x(10)*x(11) + 1/J_z*tau_psi;

%------------- END OF CODE --------------