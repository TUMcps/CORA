function dx = vanderPolEq(x,u)
% vanderPolEq - system dynamics for the Van-der-Pol oscillator 
%               (see Sec. VII in [1])
%
% Syntax:
%    dx = vanderPolEq(x,u)
%
% Inputs:
%    x - state vector
%    u - input vector
%
% Outputs:
%    dx - time-derivate of the system state
% 
% References:
%    [1] M. Althoff et al. "Reachability analysis of nonlinear systems with 
%        uncertain parameters using conservative linearization", CDC 2008

% Authors:       Matthias Althoff
% Written:       22-May-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    mu = 1;

    dx(1,1) = x(2);
    dx(2,1) = mu*(1-x(1)^2)*x(2)-x(1)+u(1);
    
% ------------------------------ END OF CODE ------------------------------
