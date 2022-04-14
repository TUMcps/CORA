function dx = lotkaVolterraCont(x,u)
% lotkaVolterraCont - Lotka-Volterra model from Example 5.2.3 in [1]
%    ... suffix Cont for continuous version (as opposed to standard hybrid)
%
% Syntax:  
%    dx = lotkaVolterraCont(x,u)
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

dx(1,1) = x(1)*(1 - (x(1) + 0.85*x(2) + 0.5*x(5)));
dx(2,1) = x(2)*(1 - (x(2) + 0.85*x(3) + 0.5*x(1)));
dx(3,1) = x(3)*(1 - (x(3) + 0.85*x(4) + 0.5*x(2)));
dx(4,1) = x(4)*(1 - (x(4) + 0.85*x(5) + 0.5*x(3)));
dx(5,1) = x(5)*(1 - (x(5) + 0.85*x(1) + 0.5*x(4)));

end

%------------- END OF CODE --------------