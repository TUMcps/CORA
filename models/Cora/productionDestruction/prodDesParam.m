function f = prodDesParam(x,u)
% prodDesParam - system dynamics for the production-destruction benchmark
%                (see Sec. 3.1 in [1])
%
% Syntax:  
%    f = prodDesParam(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
%    p - parameter vector
%
% Outputs:
%    f - time-derivate of the system state
% 
% References:
%    [1] L. Geretti, “ARCH-COMP20 Category Report: Continuous and Hybrid 
%        Systems with Nonlinear Dynamics", 2020

% Author:        Matthias Althoff
% Written:       22-May-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

    f(1,1) = (-x(1)*x(2))/(1+x(1));
    f(2,1) = (x(1)*x(2))/(1+x(1)) - x(4) * x(2);
    f(3,1) = x(4) * x(2);
    f(4,1) = 0;

end

%------------- END OF CODE --------------