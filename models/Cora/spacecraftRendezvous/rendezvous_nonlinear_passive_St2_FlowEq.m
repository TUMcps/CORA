function dx = rendezvous_nonlinear_passive_St2_FlowEq(x,u)
% rendezvous_nonlinear_passive_St2_FlowEq - system dynamics for mode 2 for
%    the spacecraft rendezvous benchmark described in [1] 
%
% Syntax:  
%    dx = rendezvous_nonlinear_passive_St2_FlowEq(x,u)
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

dx(3,1) = 0.0002624*x(2) - 0.576038456806855922248709263104*x(1) - 19.2299796*x(3) + 0.008750586984672*x(4) - (1434960000000000000.0*(x(1) + 42164000))/((x(1) + 42164000)^2 - x(2)^2)^(3/2) + 807.15359572684597539321366928407;

dx(4,1) = - 0.0002624*x(1) - 0.575980856806855922248709263104*x(2) - 0.008750586984672*x(3) - 19.2299766*x(4) - (1434960000000000000.0*x(2))/((x(1) + 42164000)^2 - x(2)^2)^(3/2);

dx(5,1) = 1;

%------------- END OF CODE --------------