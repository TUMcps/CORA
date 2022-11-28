function f = doublePendulum(x,u)
% doublePendulum - double pendulum system from Section 3.6 in [1]
%
% Syntax:  
%    f = doublePendulum(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector (= T)
% 
% Outputs:
%    f - time-derivate of the system state
%
% Reference:
%   [1] T. Johnson and et al. "ARCH-COMP21 Category Report: Artificial 
%       Intelligence and Neural Network Control Systems (AINNCS) for 
%       Continuous and Hybrid Systems Plants", 2021

% Author:       Niklas Kochdumper
% Written:      08-December-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parameter
    m = 0.5; 
    L = 0.5; 
    c = 0;
    g = 1;

    % dynamic equations
    s = cos(x(1)-x(2));

    f = [x(3); ...
         x(4); ...
         (s*(x(3)^2*sin(x(1)-x(2)) - s*(g/L*sin(x(1)) - ...
          0.5*x(4)^2*sin(x(1)-x(2)) + (u(1)-c*x(3))/(2*L^2*m)) + ...
          g/L*sin(x(2)) + (u(2)-c*x(4))/(L^2*m))) ...
          /(2*(0.5*cos(x(1)-x(2))^2-1)) - 0.5*x(4)^2*sin(x(1)-x(2)) + ...
          g/L*sin(x(1)) + (u(1)-c*x(3))/(2*L^2*m); ...
          ...
          -(x(3)^2*sin(x(1)-x(2)) - s*(g/L*sin(x(1)) - ...
          0.5*x(4)^2*sin(x(1)-x(2)) + (u(1)-c*x(3))/(2*L^2*m)) + ...
          g/L*sin(x(2)) + (u(2)-c*x(4))/(L^2*m) )/ ...
          (0.5*cos(x(1)-x(2))^2 - 1)];
     
end

%------------- END OF CODE --------------