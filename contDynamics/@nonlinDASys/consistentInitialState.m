function y0 = consistentInitialState(obj,x0,y0,u0)
% consistentInitialState - returns a consistent initial algebraic state,
%    i.e., a value for y(0) for which all algebraic equations 0 = g(x,y,u)
%    are fulfilled
%
% Syntax:
%    y0 = consistentInitialState(obj,x0,y0,u0)
%
% Inputs:
%    x0 - initial dynamic state
%    y0 - guess for the initial algebraic state
%    u0 - initial input
%
% Outputs:
%    y0 - updated initial algebraic state
%
% Example: 
%    -

% Authors:       Matthias Althoff
% Written:       18-August-2016
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%init
converged = false;

while ~converged
    l = obj.conFile(x0, y0, u0);
    [~,~,~,~,~,F] = obj.jacobian(x0, y0, u0);

    %evaluate jacobian
    delta_y = F\(-l);
    
    %check convergence
    if norm(delta_y) < 1e-10
        converged = true;
    end
    
    y0 = y0 + delta_y;
end

% ------------------------------ END OF CODE ------------------------------
