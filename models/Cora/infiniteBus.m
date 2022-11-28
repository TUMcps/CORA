function f = infiniteBus(x,u)
% infiniteBus - infinite-bus power system from [1]
%
% Syntax:  
%    f = infiniteBus(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
%
% Outputs:
%    f - time-derivate of the system state
%
% Reference:
%    [1] M. Althoff. "Formal Verification of Power Systems", submitted to
%        ARCH 2022.

% Author:       Matthias Althoff
% Written:      02-Jun-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    f = [x(2);...
         15*pi*(1 - 5*sin(x(1)) - 0.04*x(2))];
          
end

%------------- END OF CODE --------------