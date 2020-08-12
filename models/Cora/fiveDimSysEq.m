function dx = fiveDimSysEq(x,u)
% fiveDimSysEq - linear system from [1, Sec. 3.2.3] modelled as a nonlinear
%                system (required for unit test)
%
% Syntax:  
%    dx = fiveDimSysEq(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
%
% Outputs:
%    dx - time-derivate of the system state
% 
% References:
%    [1] M. Althoff, “Reachability analysis and its application to the 
%        safety assessment of autonomous cars", Dissertation, TUM 2010

% Author:        Niklas Kochdumper
% Written:       19-June-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

    dx(1,1) = -x(1) - 4*x(2) + u(1); 
    dx(2,1) = 4*x(1) - x(2) + u(2); 
    dx(3,1) = -3*x(3) + x(4) + u(3); 
    dx(4,1) = -x(3) -3*x(4) + u(4); 
    dx(5,1) = -2*x(5) + u(5); 
    
%------------- END OF CODE --------------