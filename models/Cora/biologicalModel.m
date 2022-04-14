function dx = biologicalModel(x,u)
% biologicalModel - Biological Model 1 from Example 5.2.4 in [1]
%
% Syntax:  
%    dx = biologicalModel(x,u)
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

dx(1,1) = -0.4*x(1) - x(1)*x(6);
dx(2,1) = 0.4*x(1) - x(2);
dx(3,1) = x(2) - 5*x(3)*x(4);
dx(4,1) = 5*x(5)*x(6) - 5*x(3)*x(4);
dx(5,1) = -5*x(5)*x(6) + 5*x(3)*x(4);
dx(6,1) = 0.5*x(7) - 5*x(5)*x(6);
dx(7,1) = -0.5*x(7) + 5*x(5)*x(6);

end