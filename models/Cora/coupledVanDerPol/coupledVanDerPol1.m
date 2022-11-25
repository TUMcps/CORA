function f = coupledVanDerPol1(x,u)
% coupledVanDerPol1 - system dynamics for the coupled Van-der-Pol 
%                     oscillator with stiffness 1 
%                     (see Example 5.2.2 in [1])
%
% Syntax:  
%    f = coupledVanDerPol1(x,u)
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

% Author:        Matthias Althoff
% Written:       22-May-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

    mu = 1;
    
    f(1,1) = x(2);
    f(2,1) = mu*(1-x(1)^2)*x(2) - 2*x(1) + x(3);
    f(3,1) = x(4);
    f(4,1) = mu*(1-x(3)^2)*x(4) - 2*x(3) + x(1);
end

%------------- END OF CODE --------------