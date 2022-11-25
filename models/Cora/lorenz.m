function dx = lorenz(x,u)
% lorenz - Lorenz system from Example 3.4.2 in [1]
%    note: this system is known as the Lorenz attractor,
%          it is a well-known strange attractor 
%
% Syntax:  
%    dx = lorenz(x,u)
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
% Written:       04-February-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

s = 10;
p = 8/3;
b = 28;

dx(1,1) = s * (x(2) - x(1));
dx(2,1) = x(1) * (p - x(3)) - x(2);
dx(3,1) = x(1) * x(2) - b * x(3);
    
end

%------------- END OF CODE --------------