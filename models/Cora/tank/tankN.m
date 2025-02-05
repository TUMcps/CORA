function f = tankN(x,u,h,n_x)
% tankN - system dynamics for the discrete-time version of the tank 
%   benchmark with n_x states (see Sec. VII in [1])
%
% Syntax:
%    f = tankN(x,u,h)
%
% Inputs:
%    x - state vector
%    u - input vector
%    h - time step size
%
% Outputs:
%    f - new state vector
% 
% References:
%    [1] M. Althoff et al. "Reachability analysis of nonlinear systems with 
%        uncertain parameters using conservative linearization", CDC 2008

% Authors:       Laura Luetzow, Matthias Althoff
% Written:       22-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parameter
    k = 0.3;

    % differential equations
    f(1,1) = x(1) + h*(u(1) -    k*sqrt(x(1)));           % tank 1
    for i = 2:n_x
        if mod(i-1,3) == 0
            eval(sprintf("f(%d,1)=x(%d)+h*(u(%d)+k*(sqrt(x(%d)) - sqrt(x(%d))));", ...
                i,i,(i-1)/3+1,i-1,i));
        else
            eval(sprintf("f(%d,1)=x(%d)+h*(k*(sqrt(x(%d)) - sqrt(x(%d))));", ...
                i,i,i-1,i));
        end
    end


% ------------------------------ END OF CODE ------------------------------
