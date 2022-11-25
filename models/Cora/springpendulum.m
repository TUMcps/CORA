function dx = springpendulum(x,u)
% springpendulum - spring-pendulum model from [1, Ex. 3.3.12]
% 
%           name    | init set [1]  | meaning
%   x(1):   r       | 1.2           | length of the spring
%   x(2):   theta   | 0.5           | angle between the spring and the vertical
%   x(3):   v_r     | 0             | radial velocity
%   x(4):   v_theta | 0             | angular velocity
%
% Syntax:  
%    dx = springpendulum(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
%
% Outputs:
%    dx - time-derivate of the system state
%
% Reference:
%   [1] X. Chen. "Reachability Analysis of Non-Linear Hybrid Systems Using
%       Taylor Models"

% Author:        Mark Wetzlinger
% Written:       12-October-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

g = 9.81;   % gravity
k = 2;      % spring constant
L = 1;      % natural length of the spring

dx(1,1) = x(3);
dx(2,1) = x(4);
dx(3,1) = x(1)*x(4)^2 + g*cos(x(2)) - k*(x(1) - L);
dx(4,1) = -( 2*x(3)*x(4) + g*sin(x(2)) ) / x(1);
    
end

%------------- END OF CODE --------------