function dx = roessler(x,u)
% roessler - Roessler attractor from Example 3.4.3 in [1]
%    note: a similar model is the Lorenz attractor (lorenz.m)
%
% Syntax:  
%    dx = roessler(x,u)
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

a = 0.2;
b = 0.2;
c = 5.7;

dx(1,1) = -x(2) - x(3);
dx(2,1) = x(1) + a * x(2);
dx(3,1) = b + x(3)*(x(1) - c);
    
end

%------------- END OF CODE --------------