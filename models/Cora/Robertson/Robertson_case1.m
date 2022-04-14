function f = Robertson_case1(x,u)
% Robertson_case1 - Robertson chemical reaction benchmark system from [1]
%
% Syntax:  
%    f = Robertson_case1(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
%
% Outputs:
%    f - time-derivate of the system state
%
% Reference:
%   [1] H. H. Robertson. The solution of a set of reaction rate equations.
%       In "Numerical analysis: an introduction", page 178-182.
%       Academic Press, 1966.

% Author:       Mark Wetzlinger
% Written:      17-May-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

alpha = 0.4;
beta = 1e2;
gamma = 1e3;

f = [-alpha*x(1) + beta*x(2)*x(3);...
     alpha*x(1) - beta*x(2)*x(3) - gamma*x(2)^2;...
     gamma*x(2)^2];
     
end

%------------- END OF CODE --------------