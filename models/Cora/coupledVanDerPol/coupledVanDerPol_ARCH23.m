function dx = coupledVanDerPol_ARCH23(x,u)
% coupledVanDerPol_ARCH23 - system dynamics for the Van-der-Pol oscillator,
%    taken from the ARCH23 competition
%
% Syntax:
%    dx = coupledVanDerPol_ARCH23(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
%
% Outputs:
%    dx - time-derivate of the system state
% 

% Authors:       Mark Wetzlinger
% Written:       06-July-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% x1' == y1 &amp;
% y1' == mu*(1-x1^2)*y1 + b*(x2-x1) - x1 &amp;
% x2' == y2 &amp;
% y2' == mu*(1-x2^2)*y2 - b*(x2-x1) - x2&amp;
% b' == 0

% mu = 1;

dx = x(2);
% dx(2,1) = mu*(1-x(1)^2)*x(2) + x(5)*(x(3)-x(1)) - x(1);
dx(2,1) = (1-x(1)^2)*x(2) + x(5)*(x(3)-x(1)) - x(1);
dx(3,1) = x(4);
% dx(4,1) = mu*(1-x(3)^2)*x(4) - x(5)*(x(3)-x(1)) - x(3);
dx(4,1) = (1-x(3)^2)*x(4) - x(5)*(x(3)-x(1)) - x(3);
dx(5,1) = 0;
    
% ------------------------------ END OF CODE ------------------------------
