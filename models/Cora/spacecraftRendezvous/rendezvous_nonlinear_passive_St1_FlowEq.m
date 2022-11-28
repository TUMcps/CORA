function dx = rendezvous_nonlinear_passive_St1_FlowEq(x,u)
% rendezvous_nonlinear_passive_St1_FlowEq - system dynamics for mode 1 for
%    the spacecraft rendezvous benchmark described in [1] 
%
% Syntax:  
%    dx = rendezvous_nonlinear_passive_St1_FlowEq(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
%
% Outputs:
%    dx - time-derivate of the system state
%
% References: 
%   [1] N. Chan et al. "Verifying safety of an autonomous spacecraft 
%       rendezvous mission (Benchmark proposal)"  

% Author:       Niklas Kochdumper
% Written:      23-December-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


dx(1,1) = x(3);

dx(2,1) = x(4);

dx(3,1) = 0.000201*x(2) - 0.057638256806855922248709263104*x(1) - 2.8999508*x(3) + 0.008759786984672*x(4) - (1434960000000000000.0*(x(1) + 42164000))/((x(1) + 42164000)^2 - x(2)^2)^(3/2) + 807.15359572684597539321366928407;

dx(4,1) = - 0.000174*x(1) - 0.066493256806855922248709263104*x(2) - 0.008741346984672*x(3) - 2.9030026*x(4) - (1434960000000000000.0*x(2))/((x(1) + 42164000)^2 - x(2)^2)^(3/2);

dx(5,1) = 1;

%------------- END OF CODE --------------