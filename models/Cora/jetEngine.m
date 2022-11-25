function dx = jetEngine(x,u)
% jetEngine - jet engine model from [1], Example 7
%
% Syntax:  
%    dx = jetEngine(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
%
% Outputs:
%    dx - time-derivate of the system state
% 
% References:
%    [1] Under-approximate Flowpipes for Non-linear Continuous Systems
%     Xin Chen, Sriram Sankaranarayanan, Erika Abraham
%     ISBN: 978-0-9835678-4-4.

% Author:        Mark Wetzlinger
% Written:       31-January-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

dx(1,1) = -x(2) - 1.5*x(1)^2 - 0.5*x(1)^3 - 0.5;
dx(2,1) = 3*x(1) - x(2);

end

%------------- END OF CODE --------------