function f = dynamics_approaching(x,u)
% dynamics_approaching - system dynamics for mode "approaching" for the 
%                        spacecraft rendevous benchmark described in [1] 
%
% Syntax:  
%    f = dynamics_approaching(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
%
% Outputs:
%    f - time-derivate of the system state
%
% References: 
%   [1] N. Chan et al. "Verifying safety of an autonomous spacecraft 
%       rendezvous mission (Benchmark proposal)"  

% Author:       Niklas Kochdumper
% Written:      23-December-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

f(1,1) = x(3);

f(2,1) = x(4);

f(3,1) = 0.000201*x(2) - 0.057638256806855922248709263104*x(1) - ... 
         2.8999508*x(3) + 0.008759786984672*x(4) - ...
         (1434960000000000000.0*(x(1) + 42164000))/((x(1) + ...
         42164000)^2 - x(2)^2)^(3/2) + 807.15359572684597539321366928407;

f(4,1) = - 0.000174*x(1) - 0.066493256806855922248709263104*x(2) - ...
         0.008741346984672*x(3) - 2.9030026*x(4) - ...
         (1434960000000000000.0*x(2))/((x(1) + 42164000)^2 - x(2)^2)^(3/2);

end

%------------- END OF CODE --------------