function f = dynamics_attempt(x,u)
% dynamics_attempt - system dynamics for mode "rendezvous attempt" for the 
%                    spacecraft rendevous benchmark described in [1] 
%
% Syntax:  
%    f = dynamics_attemptx,u)
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

f(3,1) = 0.0002624*x(2) - 0.576038456806855922248709263104*x(1) - ...
         19.2299796*x(3) + 0.008750586984672*x(4) - ...
         (1434960000000000000.0*(x(1) + 42164000))/((x(1) + ...
         42164000)^2 - x(2)^2)^(3/2) + 807.15359572684597539321366928407;

f(4,1) = - 0.0002624*x(1) - 0.575980856806855922248709263104*x(2) - ...
         0.008750586984672*x(3) - 19.2299766*x(4) - ...
         (1434960000000000000.0*x(2))/((x(1) + 42164000)^2 - x(2)^2)^(3/2);

end

%------------- END OF CODE --------------