function y0 = consistentInitialState(obj, x0, y0, u0)
% consistentInitialState - returns a consisten initial state, i.e. one for
% which the algebraic equations are fulfilled
%
% Syntax:  
%    example_nonlinearDA_reach_01_powerSystem_3bus()
%
% Inputs:
%    x0 - initial dynamic state
%    y0 - guessed initial algeraic state (changed by this function)
%    u0 - initial input
%
% Outputs:
%    y0 - updated initial algebraic state
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      18-August-2016
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

%init
converged = 0;

while ~converged
    l = obj.conFile(x0, y0, u0);
    [~,~,~,~,~,F] = obj.jacobian(x0, y0, u0);

    %evaluate jacobian
    delta_y = F\(-l);
    
    %check convergence
    if norm(delta_y)<1e-10
        converged = 1;
    end
    
    y0 = y0 + delta_y;
end

%------------- END OF CODE --------------
        