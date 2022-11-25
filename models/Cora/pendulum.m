function dx = pendulum(x,u,P)
% laubLoomis - dynamic equation for the a swinging pendulum
%
% Syntax:  
%    dx = pendulum(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
%    P - parameters (length of the rod and friction coefficient)
%
% Outputs:
%    dx - time-derivative of the system state
% 
% References:
%    -

% Author:        Mark Wetzlinger
% Written:       03-August-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% constants
g = 9.81;       % earth acceleration (m/s^2)

% differential equations
dx(1,1) = x(2);
dx(2,1) = -P.mu*x(2) - g/P.L * sin(x(1));

%------------- END OF CODE --------------