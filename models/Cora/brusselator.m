function f = brusselator(x,u)
% brusselator - brusselator system from Example 3.4.1 in [1]
%
% Syntax:  
%    f = brusselator(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
%
% Outputs:
%    f - time-derivate of the system state
%
% Reference:
%   [1] X. Chen. "Reachability Analysis of Non-Linear Hybrid Systems Using
%       Taylor Models"

% Author:       Niklas Kochdumper
% Written:      19-June-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    f = [1 + x(1)^2*x(2) - 2.5*x(1);...
         1.5*x(1) - x(1)^2*x(2)];
     
end

%------------- END OF CODE --------------