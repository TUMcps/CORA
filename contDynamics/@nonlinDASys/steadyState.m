function [x0,y0] = steadyState(obj, x0, y0, u0)
% steadyState - returns a steady state of a DAE system based on a
%    Newton-Raphson iteration if one exists
%
% Syntax:
%    [x0,y0] = steadyState(obj, x0, y0, u0)
%
% Inputs:
%    x0 - guessed steady state of dynamic variables
%    y0 - guessed steady state of algebraic variables
%    u0 - input
%
% Outputs:
%    x0 - steady state of dynamic variables
%    y0 - steady state of algebraic variables
%
% Example: 
%    -

% Authors:       Matthias Althoff
% Written:       08-June-2022
% Last update:   23-November-2022 (MW, call derivatives for Jacobian)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if jacobian has been computed, otherwise call derivatives
try
    nargin(obj.jacobian);
catch
    opt.tensorOrder = 1;
    derivatives(obj,opt);
end

%init
converged = false;

while ~converged
    % returned value of dynamic equations
    k = obj.dynFile(x0, y0, u0);
    % returned value of algebraic equations
    l = obj.conFile(x0, y0, u0);
    % Jacobians
    [A,~,C,D,~,F] = obj.jacobian(x0, y0, u0);

    %evaluate Jacobian
    M = [A C; D F];
    delta_z = M\(-[k;l]);
    
    %check convergence
    if norm(delta_z) < 1e-10
        converged = true;
    end
    
    %update steady state solution
    %split into delta_x and delta_y
    delta_x = delta_z(1 : length(x0));
    delta_y = delta_z((length(x0)+1) : end);
    
    %update x0
    x0 = x0 + delta_x;
    %update y0
    y0 = y0 + delta_y;
end

% ------------------------------ END OF CODE ------------------------------
